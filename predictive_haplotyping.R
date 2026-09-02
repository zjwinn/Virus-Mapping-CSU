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
  dplyr::filter(trait == "Virus") |>
  dplyr::filter(break_bft == 1) |>
  dplyr::filter(chr_wheat != "3B")

gwas_results_predcs <- read.csv(
  file = file.path(getwd(), "PreDCS_GWAS/gwas_results_continuous_pre_yd.csv")
) |>
  dplyr::filter(trait == "WSMV") |>
  dplyr::filter(break_bft == 1)
gwas_results <- dplyr::bind_rows(
  gwas_results,
  gwas_results_predcs
)
# Get top pos per chr_wheat
subset_by_top_chr_pos <- gwas_results |>
  dplyr::group_by(chr_wheat) |>
  dplyr::slice_max(order_by = neg_log_pval, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(chr_wheat, snp_id, pos) |>
  dplyr::semi_join(gwas_results, by = c("chr_wheat", "pos"))

# Read in genetic data
geno <- gaston::read.vcf(
  file.path(getwd(), "virus_mapping_2026_production_final.vcf.gz"),
  convert.chr = FALSE
)

# Subset to only post
postdcs_blues <- read.csv(
  file.path(getwd(), "multiple_environment_adjustment/blues_multiple_environment_post_yd.csv")
)$rID |>
  unique()

# Pull genotypic information and look at frequency
geno_sub <- geno[
  geno@ped$id %in% postdcs_blues,
  geno@snps$id %in% c("S3B_23933338", "S6D_6310865")
]
geno_sub_mat <- gaston::as.matrix(geno_sub)
colnames(geno_sub_mat) <- geno_sub@snps$id

# Genotype frequencies
round(table(geno_sub_mat[, 1]) / nrow(geno_sub), 2)
round(table(geno_sub_mat[, 2]) / nrow(geno_sub), 2)

# Minor allele frequencies (MAF)
p_af <- colMeans(geno_sub_mat, na.rm = TRUE) / 2
p_maf <- pmin(p_af, 1 - p_af)
round(p_maf, 2)

# Reassign chromosome names
geno@snps$chr <- gsub("Chr", "", geno@snps$chr)

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

  if (i == 3) {
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
  best_k <- as.character(sil_analysis$data$clusters[which.max(sil_analysis$data$y)]) |> as.numeric()

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

# Just pull 1mbp in either direction for 7D to address erratic LD plotting
mat_7D <- gaston::as.matrix(
  geno[
    ,
    (geno@snps$chr == "7D" &
      geno@snps$pos >= 5000000 &
      geno@snps$pos <= 10000000)
  ]
)

# Pull region for Wsm2 and Cmc4
mat_3B <- gaston::as.matrix(geno[, geno@snps$id %in% ld_per_chr$`3B`$Var2])
mat_6D <- gaston::as.matrix(geno[, geno@snps$id %in% ld_per_chr$`6D`$Var2])

# Execute the pipeline
set.seed(529)
results_3B <- run_silhouette_pipeline(mat_3B)
set.seed(828)
results_6D <- run_silhouette_pipeline(mat_6D)
set.seed(222)
results_4A <- run_silhouette_pipeline(mat_4A)
set.seed(915)
results_7D <- run_silhouette_pipeline(mat_7D)

# Save images
ggplot2::ggsave(
  filename = file.path(list_path, "PCA_of_Wsm2.jpg"),
  plot = results_3B$p_pca + ggplot2::labs(title = NULL),
  width = 6,
  height = 5,
  dpi = 300,
  units = "in"
)
ggplot2::ggsave(
  filename = file.path(list_path, "Sil_plot_of_Wsm2.jpg"),
  plot = results_3B$p_avg_sil + ggplot2::labs(title = NULL),
  width = 5,
  height = 5,
  dpi = 300,
  units = "in"
)
ggplot2::ggsave(
  filename = file.path(list_path, "PCA_of_Cmc4.jpg"),
  plot = results_6D$p_pca + ggplot2::labs(title = NULL),
  width = 6,
  height = 5,
  dpi = 300,
  units = "in"
)
ggplot2::ggsave(
  filename = file.path(list_path, "Sil_plot_of_Cmc4.jpg"),
  plot = results_6D$p_avg_sil + ggplot2::labs(title = NULL),
  width = 5,
  height = 5,
  dpi = 300,
  units = "in"
)
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
blues_me_postdcs <- read.csv(
  file.path(
    getwd(),
    "multiple_environment_adjustment/blues_multiple_environment_post_yd.csv"
  )
) |>
  dplyr::mutate(DATASET = "PostDCS")

blues_me_predcs <- read.csv(
  file.path(
    getwd(),
    "multiple_environment_adjustment/blues_multiple_environment_pre_yd.csv"
  )
) |>
  dplyr::mutate(DATASET = "PreDCS")

blues_me <- dplyr::bind_rows(blues_me_predcs, blues_me_postdcs) |>
  dplyr::select(
    DATASET,
    rID,
    dplyr::everything()
  )

# Extract and format cluster assignments from the previous pipeline results
qtl_3B_df <- data.frame(
  rID = results_3B$sil_df$sample,
  Wsm2 = results_3B$sil_df$cluster
)

qtl_6D_df <- data.frame(
  rID = results_6D$sil_df$sample,
  Cmc4 = results_6D$sil_df$cluster
)

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
  dplyr::left_join(qtl_3B_df, by = "rID") |>
  dplyr::left_join(qtl_6D_df, by = "rID") |>
  tidyr::drop_na(qtl_4A, qtl_7D, Wsm2, Cmc4)

# Calculate GRM
grm <- sommer::A.mat(
  gaston::as.matrix(geno[geno@ped$id %in% blues_me$rID, ]) - 1,
  min.MAF = 0.05
)

# Calculate PCA
pcs <- prcomp(
  gaston::as.matrix(geno[geno@ped$id %in% blues_me$rID, ]) - 1
)

# Extract PCs
pcs_dataset <- pcs$x |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "rID") |>
  as.data.frame()

# Now add markers to full dataframe
mat_3B_df <- mat_3B |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "rID") |>
  as.data.frame()
mat_4A_df <- mat_4A |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "rID") |>
  as.data.frame()
mat_6D_df <- mat_6D |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "rID") |>
  as.data.frame()
mat_7D_df <- mat_7D |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "rID") |>
  as.data.frame()

# Now add to blues_Me
blues_me_haplo <- blues_me |>
  dplyr::left_join(pcs_dataset[, 1:11], by = "rID") |>
  dplyr::left_join(mat_3B_df, by = "rID") |>
  dplyr::left_join(mat_4A_df, by = "rID") |>
  dplyr::left_join(mat_6D_df, by = "rID") |>
  dplyr::left_join(mat_7D_df, by = "rID")

# Initialize a list to save your model fits
model_fits <- list()

# For each qtl
for (i in c("Wsm2", "qtl_4A", "Cmc4", "qtl_7D")) {
  # If it is wsm2 or cmc4
  if (i %in% c("Wsm2", "Cmc4")) {
    # Pull data and coerce types explicitly
    temp_data <- blues_me_haplo |>
      dplyr::select(
        rID,
        WSMV_predicted.value,
        WSMV_std.error,
        PC1,
        PC2,
        PC3,
        dplyr::all_of(i)
      ) |>
      tidyr::drop_na() |>
      dplyr::rename(
        y = 2,
        se = 3,
        qtl = 7
      ) |>
      dplyr::filter(rID %in% rownames(grm)) |>
      dplyr::mutate(
        rID = as.factor(rID),
        y   = as.numeric(y),
        se  = as.numeric(se),
        PC1 = as.numeric(PC1),
        PC2 = as.numeric(PC2),
        PC3 = as.numeric(PC3),
        qtl = as.factor(qtl)
      ) |>
      as.data.frame()

    # Subset GRM and force as numeric matrix
    temp_grm <- as.matrix(grm[
      rownames(grm) %in% temp_data$rID,
      colnames(grm) %in% temp_data$rID
    ])

    # Ensure GRM dimensions and order match temp_data$rID levels
    temp_grm <- temp_grm[as.character(temp_data$rID), as.character(temp_data$rID)]

    # Now run sommer mmes model
    temp_fit <- sommer::mmes(
      # fixed = y ~ PC1 + PC2 + PC3 + qtl,
      fixed = y ~ qtl,
      random = ~ sommer::vsm(sommer::ism(rID), Gu = temp_grm),
      rcov = ~units,
      W = diag(1 / (temp_data$se^2)),
      data = temp_data,
      dateWarning = FALSE,
      henderson = FALSE
    )

    # Save the fit
    model_fits[[i]] <- temp_fit
  } else if (i %in% c("qtl_4A", "qtl_7D")) {
    # Pull data
    temp_data <- blues_me_haplo |>
      dplyr::select(
        rID,
        Virus_predicted.value,
        Virus_std.error,
        PC1,
        PC2,
        PC3,
        dplyr::all_of(i)
      ) |>
      tidyr::drop_na() |>
      dplyr::rename(
        y = 2,
        se = 3,
        qtl = 7
      )

    # Subset GRM
    temp_grm <- as.matrix(grm[
      rownames(grm) %in% temp_data$rID,
      colnames(grm) %in% temp_data$rID
    ])

    # Drop rID and format
    temp_data <- temp_data |>
      dplyr::filter(rID %in% rownames(grm)) |>
      dplyr::mutate(
        rID = as.factor(rID),
        y   = as.numeric(y),
        se  = as.numeric(se),
        PC1 = as.numeric(PC1),
        PC2 = as.numeric(PC2),
        PC3 = as.numeric(PC3),
        qtl = as.factor(qtl)
      ) |>
      as.data.frame()

    # Ensure GRM dimensions and order match temp_data$rID levels
    temp_grm <- temp_grm[as.character(temp_data$rID), as.character(temp_data$rID)]

    # Now run sommer model
    temp_fit <- sommer::mmes(
      # fixed = y ~ PC1 + PC2 + PC3 + qtl,
      fixed = y ~ qtl,
      random = ~ sommer::vsm(sommer::ism(rID), Gu = temp_grm),
      rcov = ~units,
      W = diag(1 / (temp_data$se^2)),
      data = temp_data,
      dateWarning = FALSE,
      henderson = FALSE
    )

    # Save the fit
    model_fits[[i]] <- temp_fit
  } else {
    stop("QTL not recognized. Check vector i.")
  }
}

# Process model predictions, visualize, and assign coding
qtl_predictions <- list()
qtl_plots <- list()
allele_coding <- list()

# Loop
for (qtl_name in names(model_fits)) {
  temp_fit <- model_fits[[qtl_name]]

  # Predict based on the "qtl" term
  pred <- sommer::predict.mmes(temp_fit, D = "qtl")$pvals

  # Calculate CI
  pred$lower <- pred$predicted.value - qnorm(0.975) * pred$std.error
  pred$upper <- pred$predicted.value + qnorm(0.975) * pred$std.error

  # Determine R (lower Y) and S (higher Y)
  pred <- pred |> dplyr::arrange(predicted.value)

  if (nrow(pred) == 2) {
    pred$Allele <- c("R", "S")
  } else {
    # Initialize Allele vector
    pred$Allele <- NA
    pred$Allele[1] <- "R"
    pred$Allele[nrow(pred)] <- "S"

    # Use significance (approx z-test) to assign any intermediate clusters
    if (nrow(pred) > 2) {
      for (j in 2:(nrow(pred) - 1)) {
        # Test against the lowest (most Resistant)
        z_R <- (pred$predicted.value[j] - pred$predicted.value[1]) / sqrt(pred$std.error[j]^2 + pred$std.error[1]^2)
        p_R <- 2 * pnorm(-abs(z_R))

        # Test against the highest (most Susceptible)
        z_S <- (pred$predicted.value[nrow(pred)] - pred$predicted.value[j]) / sqrt(pred$std.error[nrow(pred)]^2 + pred$std.error[j]^2)
        p_S <- 2 * pnorm(-abs(z_S))

        if (p_R > 0.05 && p_S <= 0.05) {
          pred$Allele[j] <- "R" # Not sig diff from R, but sig diff from S
        } else if (p_S > 0.05 && p_R <= 0.05) {
          pred$Allele[j] <- "S" # Not sig diff from S, but sig diff from R
        } else {
          # If significantly different from both (or neither), just use the closest mean
          dist_to_R <- pred$predicted.value[j] - pred$predicted.value[1]
          dist_to_S <- pred$predicted.value[nrow(pred)] - pred$predicted.value[j]
          pred$Allele[j] <- ifelse(dist_to_R < dist_to_S, "R", "S")
        }
      }
    }
  }

  # Save predictions
  pred$QTL <- qtl_name
  qtl_predictions[[qtl_name]] <- pred

  # Plot
  p <- ggplot2::ggplot(
    pred,
    ggplot2::aes(x = as.factor(qtl), y = predicted.value)
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
      title = paste("Effect of", qtl_name),
      x = paste(qtl_name, "K-means Cluster"),
      y = "Predicted Value"
    ) +
    ggplot2::theme_classic(base_size = 13)

  qtl_plots[[qtl_name]] <- p

  # Save plot
  ggplot2::ggsave(
    filename = file.path(list_path, paste0(qtl_name, "_effect_plot.jpg")),
    plot = p,
    width = 6, height = 5, dpi = 300, units = "in"
  )

  # Store the coding mapping
  mapping <- data.frame(
    Original_Cluster = as.factor(pred$qtl),
    Coded_Allele = as.character(pred$Allele)
  )
  colnames(mapping) <- c(qtl_name, paste0(qtl_name, "_coded"))
  allele_coding[[qtl_name]] <- mapping
}

# Make vector
blues_me_haplo_final <- blues_me_haplo

# Assign the R/S coding back to datasets and relocate it next to the original QTL column
for (qtl_name in names(allele_coding)) {
  coded_col <- paste0(qtl_name, "_coded")

  blues_me_haplo_final <- blues_me_haplo_final |>
    dplyr::left_join(allele_coding[[qtl_name]], by = qtl_name) |>
    dplyr::relocate(dplyr::all_of(coded_col), .after = dplyr::all_of(qtl_name))
}

# Write out
write.csv(
  blues_me_haplo_final,
  file.path(list_path, "predictive_haplotypes_and_pheno.csv"),
  row.names = FALSE
)

# Show the Wsm2 Haplotype vs PC1, PC2, and PC3
plot_data <- blues_me_haplo_final |> tidyr::drop_na(Wsm2_coded, PC1, PC2, PC3)

p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = PC1, y = PC2, color = Wsm2_coded, shape = Wsm2_coded)) +
  ggplot2::geom_point(alpha = 0.8, size = 2) +
  ggplot2::scale_color_manual(values = c("R" = "#0072B2", "S" = "#D55E00"), name = "Wsm2 Allele") +
  ggplot2::scale_shape_manual(values = c("R" = 16, "S" = 17), name = "Wsm2 Allele") +
  ggplot2::labs(x = "PC1", y = "PC2") +
  ggplot2::theme_minimal()

p2 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = PC1, y = PC3, color = Wsm2_coded, shape = Wsm2_coded)) +
  ggplot2::geom_point(alpha = 0.8, size = 2) +
  ggplot2::scale_color_manual(values = c("R" = "#0072B2", "S" = "#D55E00"), name = "Wsm2 Allele") +
  ggplot2::scale_shape_manual(values = c("R" = 16, "S" = 17), name = "Wsm2 Allele") +
  ggplot2::labs(
    title = "Wsm2 Haplotypes vs Genome-Wide Structure",
    x = "PC1", y = "PC3"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

p3 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = PC2, y = PC3, color = Wsm2_coded, shape = Wsm2_coded)) +
  ggplot2::geom_point(alpha = 0.8, size = 2) +
  ggplot2::scale_color_manual(values = c("R" = "#0072B2", "S" = "#D55E00"), name = "Wsm2 Allele") +
  ggplot2::scale_shape_manual(values = c("R" = 16, "S" = 17), name = "Wsm2 Allele") +
  ggplot2::labs(x = "PC2", y = "PC3") +
  ggplot2::theme_minimal()

# Combine using patchwork without the + operator
p_wsm2_pca_combined <- patchwork::wrap_plots(
  p1, p2, p3,
  ncol = 3,
  guides = "collect"
)

print(p_wsm2_pca_combined)

ggplot2::ggsave(
  filename = file.path(list_path, "Wsm2_vs_GenomeWide_PCA_Multi.jpg"),
  plot = p_wsm2_pca_combined,
  width = 12, height = 4, dpi = 300, units = "in"
)

# Evaluate the interaction between Wsm2 and 4A
for (i in c("R", "S")) {
  # Pull data
  temp_data <- blues_me_haplo_final |>
    dplyr::filter(Wsm2_coded == i) |>
    tidyr::drop_na(Virus_predicted.value) |>
    dplyr::select(
      rID,
      Virus_predicted.value,
      Virus_std.error,
      qtl_4A_coded
    ) |>
    dplyr::rename(
      y = 2,
      se = 3,
      qtl = 4
    )

  # Pull GRM
  temp_grm <- grm[
    rownames(grm) %in% temp_data$rID,
    colnames(grm) %in% temp_data$rID
  ]

  # Filter data and then format columns
  temp_data <- temp_data |>
    dplyr::filter(rID %in% rownames(temp_grm)) |>
    dplyr::mutate(
      rID = as.factor(rID),
      y = as.numeric(y),
      se = as.numeric(se),
      qtl = as.factor(qtl)
    )

  # Now run sommer model
  temp_fit <- sommer::mmes(
    # fixed = y ~ PC1 + PC2 + PC3 + qtl,
    fixed = y ~ qtl,
    random = ~ sommer::vsm(sommer::ism(rID), Gu = temp_grm),
    rcov = ~units,
    W = diag(1 / (temp_data$se^2)),
    data = temp_data,
    dateWarning = FALSE,
    henderson = FALSE
  )

  # Pull predictions
  temp_pred <- sommer::predict.mmes(
    temp_fit,
    D = "qtl"
  )$pvals

  # Calculate CI
  temp_pred$lower <- temp_pred$predicted.value - qnorm(0.975) * temp_pred$std.error
  temp_pred$upper <- temp_pred$predicted.value + qnorm(0.975) * temp_pred$std.error
  temp_pred$qtl <- ifelse(temp_pred$qtl == "R", "Resistant (K=1)", "Susceptible (K=2)")

  # Plot
  p <- ggplot2::ggplot(
    temp_pred,
    ggplot2::aes(x = as.factor(qtl), y = predicted.value)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper),
      width = 0.15,
      linewidth = 0.7
    ) +
    ggplot2::geom_text(
      ggplot2::aes(y = upper, label = sprintf("%.2f", upper)),
      vjust = -0.6,
      size = 3.5
    ) +
    ggplot2::geom_text(
      ggplot2::aes(y = lower, label = sprintf("%.2f", lower)),
      vjust = 1.6,
      size = 3.5
    ) +
    ggplot2::geom_label(
      ggplot2::aes(label = sprintf("%.2f", predicted.value)),
      fill = "white",
      color = "black",
      size = 3.5,
      linewidth = 0.3
    ) +
    ggplot2::scale_y_continuous(
      limits = c(min(temp_pred$lower) - 0.1, max(temp_pred$upper) + 0.1)
    ) +
    ggplot2::labs(
      title = "Estimated Group Means of 4AL Haplotypes",
      subtitle = ifelse(
        i == "R",
        paste0("Wsm2 Positive Individuals Only (n = ", nrow(temp_data), ")"),
        paste0("Wsm2 Negative Individuals Only (n = ", nrow(temp_data), ")")
      ),
      x = "Haplotype Group",
      y = "Virus Visual Ratings (Transformed)"
    ) +
    ggplot2::theme_classic(base_size = 13)

  # Save plot
  ggplot2::ggsave(
    filename = file.path(list_path, paste0("qtl_4A_effect_plot_wsm2_", i, ".jpg")),
    plot = p,
    width = 6, height = 5, dpi = 300, units = "in"
  )
}

# Get keyfile
keyfile_path <- file.path(
  getwd(),
  "../Raw_Data_Round_2/virus_mapping_2026.tsv"
)

# Check
if (file.exists(keyfile_path)) {
  # Read in true names
  keyfile <- read.table(
    keyfile_path,
    header = TRUE
  )

  # Get names refrenced
  geno_rID_ref <- data.frame(
    rID = geno@ped$id
  ) |>
    dplyr::left_join(
      keyfile |>
        dplyr::mutate(rID = FullSampleName) |>
        dplyr::distinct(RealFullSampleName, rID),
      by = "rID"
    )

  # Pull matrix and rename
  geno_alt_mat <- geno |>
    gaston::as.matrix() |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "rID") |>
    as.data.frame() |>
    dplyr::left_join(geno_rID_ref, by = "rID") |>
    dplyr::select(rID, RealFullSampleName, dplyr::everything()) |>
    dplyr::select(-rID) |>
    tibble::column_to_rownames(var = "RealFullSampleName") |>
    as.matrix()

  # Per-marker IBS = 1 - |dosage difference| / 2, averaged over called markers
  d <- abs(geno_alt_mat["Byrd", ] - geno_alt_mat["KivariAX", ])
  ibs <- mean(1 - d / 2, na.rm = TRUE)
  print(paste("Kivari AX vs Byrd IBS =", round(ibs, 3)))

  # Read in the location information
  location_infromation <- readxl::read_excel(
    file.path(getwd(), "../Raw_Data_Round_2/Training Panel List Set 21 working.xlsx"),
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
  blues_me_haplo_final_with_ids <- blues_me_haplo_final |>
    dplyr::filter(DATASET == "PostDCS") |>
    dplyr::left_join(
      keyfile |>
        dplyr::rename(rID = FullSampleName),
      by = "rID"
    ) |>
    dplyr::select(
      1:2,
      RealFullSampleName,
      LibraryPrepId,
      dplyr::everything()
    )

  # Filter for the double resistant and double susceptible haplotypes
  blues_me_filter <- blues_me_haplo_final_with_ids |>
    dplyr::filter(
      (qtl_4A_coded == "R" & qtl_7D_coded == "R") |
        (qtl_4A_coded == "S" & qtl_7D_coded == "S")
    ) |>
    dplyr::filter(startsWith(as.character(LibraryPrepId), "2025")) |>
    dplyr::filter(RealFullSampleName %in% location_infromation$FullSampleName) |>
    dplyr::distinct(RealFullSampleName, .keep_all = TRUE) |>
    dplyr::arrange(Virus_predicted.value)

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
  res_check <- "BrawlCLPlus"

  # Pull checks
  blues_me_checks <- blues_me_filter |>
    dplyr::filter(RealFullSampleName == sus_check | RealFullSampleName == res_check) |>
    dplyr::bind_rows(
      data.frame(
        RealFullSampleName = c("Fortress+BrawlCLPlus", "NTC")
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
    dplyr::select(
      367:383,
      dplyr::everything()
    )

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
    # Pull data
    temp_data <- blues_me_haplo_final |>
      dplyr::filter(DATASET == "PostDCS") |>
      tidyr::drop_na(Virus_predicted.value) |>
      dplyr::select(
        rID,
        Virus_predicted.value,
        Virus_std.error,
        qtl_4A_coded
      ) |>
      dplyr::rename(
        y = 2,
        se = 3,
        qtl = 4
      )

    # Pull GRM
    temp_grm <- grm[
      rownames(grm) %in% temp_data$rID,
      colnames(grm) %in% temp_data$rID
    ]

    # Filter data and then format columns
    temp_data <- temp_data |>
      dplyr::filter(rID %in% rownames(temp_grm)) |>
      dplyr::mutate(
        rID = as.factor(rID),
        y = as.numeric(y),
        se = as.numeric(se),
        qtl = as.factor(qtl)
      )

    # Fit model
    temp_fit <- sommer::mmes(
      # fixed = y ~ PC1 + PC2 + PC3 + qtl,
      fixed = y ~ qtl,
      random = ~ sommer::vsm(sommer::ism(rID), Gu = temp_grm),
      rcov = ~units,
      W = diag(1 / (temp_data$se^2)),
      data = temp_data,
      dateWarning = FALSE,
      henderson = FALSE
    )

    # Save RDS
    saveRDS(temp_fit, file.path(list_path, "4a_effect_model.RDS"))
  }

  # Now load
  temp_fit <- readRDS(temp_check)

  # Now predict
  temp_pred <- sommer::predict.mmes(temp_fit, D = "qtl")$pvals

  # Calculate CI
  temp_pred$lower <- temp_pred$predicted.value - qnorm(0.975) * temp_pred$std.error
  temp_pred$upper <- temp_pred$predicted.value + qnorm(0.975) * temp_pred$std.error

  # Recode group labels
  temp_pred$qtl <- factor(
    temp_pred$qtl,
    levels = c("R", "S"),
    labels = c("Resistant (K=1)", "Susceptible (K=2)")
  )

  # Plot
  p <- ggplot2::ggplot(
    temp_pred,
    ggplot2::aes(x = qtl, y = predicted.value)
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

  # Now get real names and their allelic state
  table_for_pub <- blues_me_haplo_final_with_ids |>
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
      !startsWith(RealFullSampleName, "AP") &
      !startsWith(RealFullSampleName, "BASF") &
      !startsWith(RealFullSampleName, "LCH") &
      !startsWith(RealFullSampleName, "TAM") &
      !grepl("^[0-9]", RealFullSampleName)) |>
    dplyr::distinct() |>
    dplyr::mutate(
      RealFullSampleName = gsub("CL", " CL", RealFullSampleName),
      RealFullSampleName = gsub("AX", " AX", RealFullSampleName),
      RealFullSampleName = gsub("SF", " SF", RealFullSampleName),
      RealFullSampleName = gsub("Plus", " Plus", RealFullSampleName),
      RealFullSampleName = gsub("2.0", " 2.0", RealFullSampleName)
    ) |>
    dplyr::rename(
      `Cultivar Name` = RealFullSampleName,
      `K-means Cluster` = qtl_4A,
      `Predicted Allelic State` = qtl_4A_coded,
      `Visual Virus Rating (Transformed)` = Virus_predicted.value,
      `SE` = Virus_std.error
    ) |>
    dplyr::arrange(`Predicted Allelic State`, `Visual Virus Rating (Transformed)`) |>
    dplyr::filter(
      !`Cultivar Name` %in% c("Fortify SF", "Snowmass 2.0")
    )

  # Write out table
  write.csv(
    table_for_pub,
    file.path(list_path, "cultivar_haps.csv"),
    row.names = FALSE
  )
}
