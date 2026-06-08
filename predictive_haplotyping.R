# Create folder for list path
list_path <- file.path(getwd(), "predictive_haplotyping")
dir.create(list_path, showWarnings = FALSE)

# Set seed
set.seed(222)

# Okabe-Ito colorblind-safe palette
cb_palette <- c(
  "#D55E00", "#56B4E9", "#009E73", "#E69F00",
  "#0072B2", "#CC79A7", "#F0E442", "#000000"
)

# Read and keep break_bft == 1
gwas_results <- read.csv(
  file.path(getwd(), "PostDCS_GWAS/gwas_results_continuous_post_yd.csv")
) |>
  dplyr::filter(break_bft == 1)

# Get top pos per chr_wheat
subset_by_top_chr_pos <- gwas_results |>
  dplyr::group_by(chr_wheat) |>
  dplyr::slice_max(order_by = neg_log_pval, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(chr_wheat, snp_id, pos) |>
  dplyr::semi_join(gwas_results, by = c("chr_wheat", "pos"))

# Read in genetic data
geno <- gaston::read.vcf(
  file.path(getwd(), "virus_gwas_2025_production_final.vcf.gz"),
  convert.chr = FALSE
)

# LD threshold for labeling
cutoff <- 0.2

# Get the length of LD distance per chromosome
ld_per_chr <- list()

# Loop over top hits and plot LD for the focal SNP
for (i in seq_len(nrow(subset_by_top_chr_pos))) {
  # Pull data for this top hit
  temp_data <- subset_by_top_chr_pos[i, , drop = FALSE]
  focal_snp <- as.character(temp_data$snp_id)
  focal_chr <- as.character(temp_data$chr_wheat)

  # Subset genotypic data to the chromosome of interest
  temp_geno <- geno[, geno@snps$chr == focal_chr]

  # Compute LD (r2) and reshape
  ld_mat <- gaston::LD(
    temp_geno,
    lim = c(1, ncol(temp_geno)),
    measure = "r2"
  )

  ld_long <- reshape2::melt(ld_mat) |>
    as.data.frame() |>
    dplyr::mutate(
      Var1 = as.character(Var1),
      Var2 = as.character(Var2)
    ) |>
    dplyr::filter(Var1 == focal_snp) |>
    dplyr::mutate(
      Position = as.numeric(sub("^[^_]*_", "", Var2)),
      r2 = value
    ) |>
    dplyr::arrange(Position)

  # Skip if no LD rows found for focal SNP
  if (nrow(ld_long) == 0) {
    message("No LD entries found for focal SNP: ", focal_snp, " (chr ", focal_chr, ")")
    next
  }

  # Ensure Position numeric and flag hits
  ld_long <- ld_long |>
    dplyr::mutate(
      Position = as.numeric(Position),
      hit = r2 >= 0.2,
      Chr = sub("_.*$", "", Var2)
    )

  # Make label for cutoff
  cutoff_label <- paste0("LD = ", cutoff)

  # Make LD plot
  p <- ggplot2::ggplot(ld_long, ggplot2::aes(x = Position / 1000000, y = r2)) +
    ggplot2::geom_point(ggplot2::aes(color = hit), alpha = 0.8, size = 1.8) +
    ggplot2::geom_hline(yintercept = cutoff, linetype = "dashed", color = "#0072B2", linewidth = 0.6) +
    ggplot2::annotate(
      "text",
      x = Inf, y = cutoff,
      label = cutoff_label,
      hjust = 1.05, vjust = -0.5,
      color = "#0072B2", size = 3
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "grey70", "TRUE" = "#D55E00"),
      labels = c("FALSE" = paste0("< ", cutoff), "TRUE" = paste0(">= ", cutoff)),
      name = "LD"
    ) +
    ggplot2::labs(
      title = paste("LD (r2) with", focal_snp),
      x = "Genomic Position (Mbp)",
      y = expression(r^2)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      strip.text = ggplot2::element_text(size = 9),
      legend.position = "right"
    )

  if (i == 2) {
    p <- p +
      ggplot2::labs(
        title = NULL
      )

    qtl_4as_ld_graph <- p
  }

  # Print the image
  print(p)

  # Save image
  ggplot2::ggsave(
    filename = paste0(list_path, "/", focal_snp, "_ld_decay_graph.jpg"),
    plot = p,
    width = 10,
    height = 5,
    dpi = 300,
    units = "in"
  )

  # Put in place
  ld_per_chr[[focal_chr]] <- ld_long |>
    dplyr::filter(hit == TRUE)

  # Remove unused objects to save memory
  remove(ld_long, ld_mat, temp_data, temp_geno, focal_chr, focal_snp, p, cutoff_label)
}

# Define a haplotype pipeline
run_silhouette_pipeline <- function(chr_matrix) {
  # Utilize factoextra for automated silhouette analysis
  k_max <- min(10, nrow(chr_matrix) - 1)
  sil_analysis <- factoextra::fviz_nbclust(chr_matrix, stats::kmeans, method = "silhouette", k.max = k_max)

  # Extract best K natively from the factoextra data structure
  best_k <- as.numeric(sil_analysis$data$clusters[which.max(sil_analysis$data$y)])

  # Apply kmeans with the optimal k
  best_km <- stats::kmeans(chr_matrix, centers = best_k, nstart = 50)
  diss <- stats::dist(chr_matrix, method = "euclidean")
  sil_obj <- cluster::silhouette(best_km$cluster, diss)

  sil_df <- data.frame(
    sample = rownames(chr_matrix),
    cluster = factor(sil_obj[, "cluster"]),
    sil_width = sil_obj[, "sil_width"]
  )

  # Run PCA
  pca <- stats::prcomp(chr_matrix, center = TRUE, scale. = TRUE)

  # Compute percent variance for labeling
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
  pc1_lab <- paste0("PC1 (", round(var_expl[1] * 100, 1), "%)")
  pc2_lab <- paste0("PC2 (", round(var_expl[2] * 100, 1), "%)")

  pca_df <- data.frame(
    Sample = rownames(chr_matrix),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    cluster = factor(best_km$cluster)
  )

  # Generate plots
  p_avg_sil <- sil_analysis +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste("Best k =", best_k))

  p_pca <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = cluster)) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::scale_color_manual(values = cb_palette) +
    ggplot2::labs(
      title = paste("PCA colored by kmeans (k =", best_k, ")"),
      x = pc1_lab,
      y = pc2_lab
    ) +
    ggplot2::theme_minimal()

  print(p_avg_sil)
  print(p_pca)

  # Return objects silently for downstream use
  return(list(
    matrix = chr_matrix,
    best_k = best_k,
    sil_analysis = sil_analysis,
    best_km = best_km,
    sil_df = sil_df,
    pca = pca,
    pca_df = pca_df,
    p_pca = p_pca,
    p_avg_sil = p_avg_sil
  ))
}

# Prepare your matrices for the targeted chromosomes
qtl_4A <- ld_per_chr$`4A`
mat_4A <- gaston::as.matrix(geno[, geno@snps$id %in% qtl_4A$Var2])

# Find the Kivari AX if you can
check <- file.path(
  getwd(),
  "../../../",
  "Annual_Breeding_Pipeline",
  "GABP_2025-2026",
  "Genetic_Information",
  "Annual_GBS_Pipeline_2025_2026_Post_QC_production_final.vcf.gz"
)

# Check if file exists
if (file.exists(check)) {
  # Read in vcf
  kivariax <- gaston::read.vcf(
    check,
    convert.chr = FALSE
  )

  # Pull Kivari
  kivariax <- kivariax[kivariax@ped$id == "KivariAX", kivariax@snps$id %in% colnames(mat_4A)] |>
    gaston::as.matrix() |>
    as.data.frame()

  # Bind rows
  mat_4A <- mat_4A |>
    as.data.frame() |>
    dplyr::bind_rows(kivariax)

  # temp_names
  temp_names <- rownames(mat_4A)

  # KNN imputation
  mat_4A <- VIM::kNN(
    mat_4A,
    k = 5,
    numFun = function(x) round(mean(x)), # keep calls on the 0/1/2 dosage scale
    imp_var = FALSE
  )

  # Apply names
  rownames(mat_4A) <- temp_names

  # Make a matrix
  mat_4A <- as.matrix(mat_4A)

  # Get IBS between Byrd and Kivari AX
  m <- gaston::read.vcf(
    check,
    convert.chr = FALSE
  )

  # Byrd and Kivari AX as a 2-row dosage matrix on matched, non-missing markers
  m <- gaston::as.matrix(m[m@ped$id %in% c("Byrd", "KivariAX"), ])

  # Per-marker IBS = 1 - |dosage difference| / 2, averaged over called markers
  d <- abs(m["Byrd", ] - m["KivariAX", ])
  ibs <- mean(1 - d / 2, na.rm = TRUE)
  print(paste("Kivari AX vs Byrd IBS =", round(ibs, 3)))
}

# Just pull 1mbp in either direction for 7D to address erratic LD plotting
mat_7D <- gaston::as.matrix(
  geno[
    ,
    (geno@snps$chr == "7D" &
      geno@snps$pos >= 5000000 &
      geno@snps$pos <= 10000000)
  ]
)

# Execute the pipeline
set.seed(222)
results_4A <- run_silhouette_pipeline(mat_4A)
set.seed(915)
results_7D <- run_silhouette_pipeline(mat_7D)

# Save images
ggplot2::ggsave(
  filename = file.path(list_path, "PCA_of_4A_QTL.jpg"),
  plot = results_4A$p_pca + ggplot2::labs(title = NULL),
  width = 6,
  height = 5,
  dpi = 300,
  units = "in"
)
ggplot2::ggsave(
  filename = file.path(list_path, "Sil_plot_of_4A_QTL.jpg"),
  plot = results_4A$p_avg_sil + ggplot2::labs(title = NULL),
  width = 5,
  height = 5,
  dpi = 300,
  units = "in"
)
ggplot2::ggsave(
  filename = file.path(list_path, "PCA_of_7D_QTL.jpg"),
  plot = results_7D$p_pca,
  width = 6,
  height = 5,
  dpi = 300,
  units = "in"
)
ggplot2::ggsave(
  filename = file.path(list_path, "Sil_plot_of_7D_QTL.jpg"),
  plot = results_7D$p_avg_sil,
  width = 5,
  height = 5,
  dpi = 300,
  units = "in"
)

# Read in across-environments BLUEs
blues_me <- read.csv(
  file.path(
    getwd(),
    "multiple_environment_adjustment/blues_multiple_environment_post_yd.csv"
  )
)

# Extract and format cluster assignments from the previous pipeline results
qtl_4A_df <- data.frame(
  rID = results_4A$sil_df$sample,
  qtl_4A = results_4A$sil_df$cluster
)

qtl_7D_df <- data.frame(
  rID = results_7D$sil_df$sample,
  qtl_7D = results_7D$sil_df$cluster
)

# Merge cluster assignments with the phenotypic data
blues_me <- blues_me |>
  dplyr::left_join(qtl_4A_df, by = "rID") |>
  dplyr::left_join(qtl_7D_df, by = "rID") |>
  dplyr::mutate(
    qtl_4A = as.factor(qtl_4A),
    qtl_7D = as.factor(qtl_7D)
  ) |>
  tidyr::drop_na(qtl_4A, qtl_7D)

# Fit the full model using the raw cluster assignments
qtl_model <- stats::lm(Virus_predicted.value ~ qtl_4A * qtl_7D, data = blues_me)

# Calculate estimated marginal means for the individual loci
emm_4A <- emmeans::emmeans(qtl_model, ~qtl_4A)
emm_7D <- emmeans::emmeans(qtl_model, ~qtl_7D)

# Generate compact letter displays for main effects natively
cld_4A <- as.data.frame(multcomp::cld(emm_4A, Letters = letters))
cld_4A$.group <- trimws(cld_4A$.group)

cld_7D <- as.data.frame(multcomp::cld(emm_7D, Letters = letters))
cld_7D$.group <- trimws(cld_7D$.group)

# Plot for Chromosome 4A main effect with CLD
p_4A <- ggplot2::ggplot(cld_4A, ggplot2::aes(x = qtl_4A, y = emmean)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  ggplot2::geom_text(ggplot2::aes(y = upper.CL, label = .group), vjust = -0.5) +
  ggplot2::labs(
    x = "Chromosome 4A Cluster",
    y = "Viral Visual Ratings (Transformed)",
    title = "Main Effect: Chromosome 4A Clusters"
  ) +
  ggplot2::theme_minimal()

# Plot for Chromosome 7D main effect with CLD
p_7D <- ggplot2::ggplot(cld_7D, ggplot2::aes(x = qtl_7D, y = emmean)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  ggplot2::geom_text(ggplot2::aes(y = upper.CL, label = .group), vjust = -0.5) +
  ggplot2::labs(
    x = "Chromosome 7D Cluster",
    y = "Viral Visual Ratings (Transformed)",
    title = "Main Effect: Chromosome 7D Clusters"
  ) +
  ggplot2::theme_minimal()

print(p_4A)
print(p_7D)

# Save images
ggplot2::ggsave(
  filename = file.path(list_path, "Clustering_Main_Effect_of_4A_QTL.jpg"),
  plot = p_4A,
  width = 6,
  height = 6,
  dpi = 300,
  units = "in"
)
ggplot2::ggsave(
  filename = file.path(list_path, "Clustering_Main_Effect_of_7D_QTL.jpg"),
  plot = p_7D,
  width = 6,
  height = 6,
  dpi = 300,
  units = "in"
)

# Calculate interaction means and generate CLD
emm_interaction <- emmeans::emmeans(qtl_model, ~ qtl_4A * qtl_7D)
cld_interaction <- as.data.frame(multcomp::cld(emm_interaction, Letters = letters))
cld_interaction$.group <- trimws(cld_interaction$.group)

# Plot for the interaction effect with CLD
p_interaction <- ggplot2::ggplot(cld_interaction, ggplot2::aes(x = qtl_4A, y = emmean, color = qtl_7D, group = qtl_7D)) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.3), size = 3) +
  ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.3)) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.3,
    position = ggplot2::position_dodge(width = 0.3)
  ) +
  ggplot2::geom_text(
    ggplot2::aes(y = upper.CL, label = .group),
    position = ggplot2::position_dodge(width = 0.3),
    vjust = -0.5,
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(values = cb_palette) +
  ggplot2::labs(
    x = "Chromosome 4A Cluster",
    y = "Viral Visual Ratings (Transformed)",
    color = "7D Cluster",
    title = "Interaction Effect: Raw 4A and 7D Clusters"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "right")

print(p_interaction)

# Save image
ggplot2::ggsave(
  filename = file.path(list_path, "Clustering_of_4A_QTL_by_7D_QTL_Interaction.jpg"),
  plot = p_interaction,
  width = 6,
  height = 6,
  dpi = 300,
  units = "in"
)

# Delineate haplotypes based on the empirical EMMs observed above
blues_me <- blues_me |>
  dplyr::mutate(
    qtl_4A_coded = dplyr::case_when(
      qtl_4A %in% c("2") ~ "R",
      .default = "S"
    ),
    qtl_7D_coded = dplyr::case_when(
      qtl_7D %in% c("1", "2", "4") ~ "R",
      .default = "S"
    )
  )

# Fit the final haplotype model testing the defined R and S main effects and interaction
haplotype_model <- stats::lm(Virus_predicted.value ~ qtl_4A_coded * qtl_7D_coded, data = blues_me)
print(summary(haplotype_model))

# Calculate estimated marginal means
haplotype_emm <- emmeans::emmeans(haplotype_model, ~ qtl_4A_coded * qtl_7D_coded)

# Generate the compact letter display
haplotype_cld <- as.data.frame(multcomp::cld(haplotype_emm, Letters = letters))
haplotype_cld$.group <- trimws(haplotype_cld$.group)

# Plot the final haplotype interaction
p_haplotype <- ggplot2::ggplot(haplotype_cld, ggplot2::aes(x = qtl_4A_coded, y = emmean, color = qtl_7D_coded, group = qtl_7D_coded)) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.3), size = 3) +
  ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.3)) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.3,
    position = ggplot2::position_dodge(width = 0.3)
  ) +
  ggplot2::geom_text(
    ggplot2::aes(y = upper.CL, label = .group),
    position = ggplot2::position_dodge(width = 0.3),
    vjust = -0.5,
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(values = cb_palette) +
  ggplot2::labs(
    x = "Chromosome 4A Haplotype",
    y = "Viral Visual Ratings (Transformed)",
    color = "7D Haplotype",
    title = "Final Haplotype Interaction: 4A and 7D"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "right")

print(p_haplotype)

# Save image
ggplot2::ggsave(
  filename = file.path(list_path, "Haplotype_Clustering_of_4A_QTL_by_7D_QTL_Interaction.jpg"),
  plot = p_haplotype,
  width = 6,
  height = 6,
  dpi = 300,
  units = "in"
)

# Now add markers to full dataframe
mat_4A_df <- mat_4A |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "rID") |>
  as.data.frame()
mat_7D_df <- mat_7D |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "rID") |>
  as.data.frame()

# Now add to blues_Me
blues_me <- blues_me |>
  dplyr::left_join(mat_4A_df, by = "rID") |>
  dplyr::left_join(mat_7D_df, by = "rID")

# Write out
write.csv(
  blues_me,
  file.path(list_path, "predictive_haplotypes_and_pheno.csv"),
  row.names = FALSE
)

# Get keyfile
keyfile_path <- file.path(
  getwd(),
  "../Raw_Data/virus_mapping_2025.tsv"
)

# Check
if (file.exists(keyfile_path)) {
  # Read in true names
  keyfile <- read.table(
    keyfile_path,
    header = TRUE
  )

  # Read in the location information
  location_infromation <- readxl::read_excel(
    file.path(getwd(), "../Raw_Data/Training Panel List Set 21 working.xlsx"),
    sheet = "Master List"
  ) |>
    dplyr::rename(
      FullSampleName = FULLSAMPLENAME,
      ROW = Row,
      COLUMN = Column,
    ) |>
    dplyr::select(
      -FLOWCELL,
      -LIBRARYPREPID,
      -LANE
    )

  # Bind in keyfile to blues_me
  blues_me <- blues_me |>
    dplyr::left_join(keyfile |> dplyr::rename(rID = FullSampleName), by = "rID")

  # Filter for the double resistant and double susceptible haplotypes
  blues_me_filter <- blues_me |>
    dplyr::filter(
      (qtl_4A_coded == "R" & qtl_7D_coded == "R") |
        (qtl_4A_coded == "S" & qtl_7D_coded == "S")
    ) |>
    dplyr::filter(startsWith(as.character(LibraryPrepId), "2025")) |>
    dplyr::filter(RealFullSampleName %in% location_infromation$FullSampleName) |>
    dplyr::distinct(RealFullSampleName, .keep_all = TRUE) |>
    dplyr::arrange(Virus_predicted.value) |>
    dplyr::select(
      113:119,
      dplyr::everything()
    )

  # Extract top 14 resistant lines (strictly R/R haplotype with the lowest virus ratings)
  blues_me_top_14 <- blues_me_filter |>
    dplyr::filter(qtl_4A_coded == "R" & qtl_7D_coded == "R") |>
    dplyr::slice_min(order_by = Virus_predicted.value, n = 14, with_ties = FALSE) |>
    dplyr::mutate(Note = "Supposedly Resistant") |>
    dplyr::select(Note, dplyr::everything())

  # Extract bottom 14 susceptible lines (strictly S/S haplotype with the highest virus ratings)
  blues_me_bottom_14 <- blues_me_filter |>
    dplyr::filter(qtl_4A_coded == "S" & qtl_7D_coded == "S") |>
    dplyr::slice_max(order_by = Virus_predicted.value, n = 14, with_ties = FALSE) |>
    dplyr::mutate(Note = "Supposedly Susceptible") |>
    dplyr::select(Note, dplyr::everything())

  # Pick positive and negative checks
  sus_check <- "Fortress"
  res_check <- "Canvas"

  # Pull checks
  blues_me_checks <- blues_me_filter |>
    dplyr::filter(RealFullSampleName == sus_check | RealFullSampleName == res_check) |>
    dplyr::bind_rows(
      data.frame(
        RealFullSampleName = c("Fortress+Canvas", "NTC")
      )
    ) |>
    dplyr::mutate(Note = c("Resistant Check", "Susceptible Check", "Artifical Het", "Non-Template Control")) |>
    dplyr::select(Note, dplyr::everything())

  # Bind
  test_plate <- rbind(
    blues_me_bottom_14,
    blues_me_top_14,
    blues_me_checks
  )

  # Make test plate with locations
  test_plate_with_loc <- test_plate |>
    dplyr::rename(FullSampleName = RealFullSampleName) |>
    dplyr::left_join(
      location_infromation,
      by = c("FullSampleName", "ROW", "COLUMN")
    ) |>
    dplyr::select(137:155, dplyr::everything())

  # Create a dummy vector of 96 line names
  line_names <- rep(test_plate_with_loc$FullSampleName, 3)

  # Convert vector to an 8x12 matrix (fills by column by default)
  plate_layout <- matrix(
    line_names,
    nrow = 8,
    ncol = 12,
    byrow = FALSE
  )

  # Assign standard 96-well plate dimension names
  rownames(plate_layout) <- LETTERS[1:8] # Rows A-H
  colnames(plate_layout) <- 1:12 # Columns 1-12

  # Now read in the other plate I am working with
  check_plate <- read.csv(
    file.path(
      getwd(),
      "../Raw_Data/selected_individuals_for_test_trays.csv"
    )
  )

  # Check if there is overlap
  all(check_plate$genotype %in% test_plate_with_loc$FullSampleName)

  # Now create a workbook for the test plate
  wb <- openxlsx::createWorkbook()

  # Add pages
  openxlsx::addWorksheet(wb, "plate_layout")
  openxlsx::addWorksheet(wb, "plate_information")

  # Write first sheet
  openxlsx::writeData(
    wb,
    sheet = "plate_layout",
    x = plate_layout,
    rowNames = TRUE,
  )

  # Write second sheet
  openxlsx::writeData(
    wb,
    sheet = "plate_information",
    x = test_plate_with_loc,
  )

  # Write out workbook
  openxlsx::saveWorkbook(wb, file.path(list_path, "Test_Plates_Map_and_Information_for_4A_and_7D_QTL.xlsx"), overwrite = TRUE)

  # Check for file
  temp_check <- file.path(list_path, "4a_effect_model.RDS")

  # If it doesn't exist
  if (!file.exists(temp_check)) {
    # Get pheno information
    temp_pheno <- blues_me |>
      dplyr::select(rID, Virus_predicted.value, Virus_std.error, qtl_4A, qtl_4A_coded) |>
      tidyr::drop_na()

    # Pull genetic information
    temp_grm <- geno[
      geno@ped$id %in% temp_pheno$rID,
      !geno@snps$id %in% colnames(mat_4A)
    ] |>
      gaston::as.matrix() - 1

    # Calculate GRM
    temp_grm <- sommer::A.mat(temp_grm)

    # Fit model
    temp_fit <- sommer::mmes(
      fixed = Virus_predicted.value ~ qtl_4A_coded,
      random = ~ sommer::vsm(sommer::ism(rID), Gu = temp_grm),
      rcov = ~units,
      W = diag(temp_pheno$std.error),
      data = temp_pheno,
      dateWarning = FALSE
    )

    # Save RDS
    saveRDS(temp_fit, file.path(list_path, "4a_effect_model.RDS"))
  }

  # Now load
  temp_fit <- readRDS(temp_check)

  # Now predict
  temp_pred <- sommer::predict.mmes(temp_fit, D = "qtl_4A_coded")$pvals

  # Calculate CI
  temp_pred$lower <- temp_pred$predicted.value - qnorm(0.975) * temp_pred$std.error
  temp_pred$upper <- temp_pred$predicted.value + qnorm(0.975) * temp_pred$std.error

  # Recode group labels
  temp_pred$qtl_4A_coded <- factor(
    temp_pred$qtl_4A_coded,
    levels = c("R", "S"),
    labels = c("Resistant (K=2)", "Susceptible (K=1)")
  )

  # Plot
  p <- ggplot2::ggplot(
    temp_pred,
    ggplot2::aes(x = qtl_4A_coded, y = predicted.value)
  ) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper),
      width = 0.15, linewidth = 0.7
    ) +
    ggplot2::geom_text(ggplot2::aes(y = upper, label = sprintf("%.2f", upper)),
      vjust = -0.6, size = 3.5
    ) +
    ggplot2::geom_text(ggplot2::aes(y = lower, label = sprintf("%.2f", lower)),
      vjust = 1.6, size = 3.5
    ) +
    ggplot2::geom_label(ggplot2::aes(label = sprintf("%.2f", predicted.value)),
      fill = "white", color = "black", size = 3.5,
      linewidth = 0.3
    ) +
    ggplot2::labs(
      x = "Haplotype Group",
      y = "Virus Visual Rating (Transformed)"
    ) +
    ggplot2::coord_cartesian(ylim = c(-0.3, 0.7)) +
    ggplot2::theme_classic(base_size = 13)

  # library
  library(patchwork)

  # Assemble
  p <- qtl_4as_ld_graph +
    (results_4A$p_pca +
      ggplot2::labs(
        title = NULL,
        color = "K Cluster"
      )) +
    p +
    plot_annotation(tag_levels = "A")

  # Save image
  ggplot2::ggsave(
    file.path(list_path, "4as_composite_image_PostDCS.jpg"),
    plot = p,
    width = 15,
    height = 5,
    units = "in",
    dpi = 300
  )

  # Show min max
  min(qtl_4A$Position) / 1000000
  max(qtl_4A$Position) / 1000000
  max(qtl_4A$Position) / 1000000 - min(qtl_4A$Position) / 1000000

  # Confirmed non-CSU entries (out-of-program checks/landrace)
  non_csu <- c("Jagalene", "Kharkof", "Scout66")

  # Pull Kivari AX
  kivariax <- data.frame(
    `Cultivar Name` = "Kivari AX",
    `K-means Cluster` = factor(
      results_4A$best_km$cluster["KivariAX"],
      levels = levels(blues_me$qtl_4A)
    ),
    `Predicted Allelic State` = ifelse(results_4A$best_km$cluster["KivariAX"] == 1, "S", "R"),
    `Visual Virus Rating (Transfomred)` = NA,
    `SE` = NA,
    check.names = FALSE
  )
  # Now get real names and their allelic state
  table_for_pub <- blues_me |>
    dplyr::select(RealFullSampleName, qtl_4A, qtl_4A_coded, Virus_predicted.value, Virus_std.error) |>
    dplyr::mutate(
      RealFullSampleName = dplyr::if_else(
        RealFullSampleName == "CO200037R", "Gabriel", RealFullSampleName
      )
    ) |>
    dplyr::filter(!RealFullSampleName %in% non_csu) |>
    dplyr::filter(!startsWith(RealFullSampleName, "CO") &
      !startsWith(RealFullSampleName, "TX") &
      !startsWith(RealFullSampleName, "OK") &
      !startsWith(RealFullSampleName, "KS") &
      !startsWith(RealFullSampleName, "NE") &
      !startsWith(RealFullSampleName, "NH") &
      !startsWith(RealFullSampleName, "MT") &
      !startsWith(RealFullSampleName, "WB") &
      !startsWith(RealFullSampleName, "NI") &
      !grepl("^[0-9]", RealFullSampleName)) |>
    dplyr::distinct() |>
    dplyr::rename(
      `Cultivar Name` = RealFullSampleName,
      `K-means Cluster` = qtl_4A,
      `Predicted Allelic State` = qtl_4A_coded,
      `Visual Virus Rating (Transfomred)` = Virus_predicted.value,
      `SE` = Virus_std.error
    ) |>
    dplyr::bind_rows(kivariax) |>
    dplyr::arrange(`Predicted Allelic State`, `Visual Virus Rating (Transfomred)`)

  # Write out table
  write.csv(
    table_for_pub,
    file.path(list_path, "cultivar_haps.csv"),
    row.names = FALSE
  )
}
