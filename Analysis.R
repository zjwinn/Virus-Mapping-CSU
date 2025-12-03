## ----setup, include=FALSE-------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----single_env_adj-------------------------------------------------------------------------------------------------------
# Library
library(ggplot2)
library(patchwork)

# Check dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("GENESIS", quietly = TRUE)) {
  BiocManager::install("GENESIS")
}

# Define the inverse normal transformation with blom offset function
INT_blom <- function(x) {
  # c = 3/8 (Blom offset)
  n <- sum(!is.na(x))
  qnorm((rank(x, na.last = "keep", ties.method = "average") - 0.375) / (n + 0.25))
}

# Read in pheno data
pheno <- read.csv("privatized_data.csv") |>
  dplyr::mutate(Y = ifelse(RESPONSE == "GY", Y * (67.25), # Conversion to Kg/H
    ifelse(RESPONSE == "TW", Y * (0.4536 / 0.3524), Y) # Conversion to Kg/Hl
  ))

# Check if the analysis needs to be ran
if (!file.exists("blues_single_environment.csv")) {
  # Get adjusted means object ready
  blues_se <- c()

  # Run adjustment within environment
  for (i in unique(pheno$YEN)) {
    # Make message
    print_message <- paste("### Analyzing YEN =", i, "###")
    border <- paste(rep("#", nchar(print_message)), collapse = "")
    print(border)
    print(print_message)
    print(border)

    # Pull data
    temp1 <- pheno[pheno$YEN == i, ]

    # For j in traits
    for (j in unique(temp1$RESPONSE)) {
      # Make message
      print_message <- paste("### Trait =", j, "###")
      border <- paste(rep("#", nchar(print_message)), collapse = "")
      print(border)
      print(print_message)
      print(border)

      # Pull trait
      temp2 <- temp1[temp1$RESPONSE == j, ]

      # Order row:col
      temp2 <- temp2[order(temp2$ROW, temp2$COLUMN), ]

      # Check if col 1 is col 1 if not set to 1
      temp2$ROW <- temp2$ROW |>
        as.factor() |>
        as.numeric()
      temp2$COLUMN <- temp2$COLUMN |>
        as.factor() |>
        as.numeric()
      temp2$Y <- as.numeric(temp2$Y)

      # Drop NA
      temp2 <- temp2 |> na.omit()

      # if temp2 has no data
      if (nrow(temp2) == 0) {
        remove(temp2)
        next
      }

      # Check type
      if (unique(temp2$DATA_TYPE) == "NORMAL") {
        # Run model in sommer
        temp3 <- try(sommer::mmes(
          fixed = Y ~ rID,
          random = ~ sommer::spl2Dc(ROW, COLUMN),
          rcov = ~units,
          data = temp2
        ))

        # If it fails toss the data
        if (class(temp3) == "try-error") {
          remove(temp2, temp3)
          next
        }

        # Calculate adjusted means
        temp4 <- sommer::predict.mmes(temp3, D = "rID")$pvals

        # Bind to adjusted means
        temp4$RESPONSE <- j
        temp4$YEN <- i
        temp4$DATA_TYPE <- "BLUE"
        blues_se <- rbind(blues_se, temp4)

        # Remove
        remove(temp2, temp3, temp4)
      } else if (unique(temp2$DATA_TYPE) == "ORDINAL") {
        # Transform y
        temp2$Y_t <- INT_blom(temp2$Y)

        # Run model in sommer
        temp3 <- try(sommer::mmes(
          fixed = Y_t ~ rID,
          random = ~ sommer::spl2Dc(ROW, COLUMN),
          rcov = ~units,
          data = temp2
        ))

        # If it fails toss the data
        if (class(temp3) == "try-error") {
          remove(temp2, temp3)
          next
        }

        # Calculate adjusted means
        temp4 <- sommer::predict.mmes(temp3, D = "rID")$pvals

        # Bind to adjusted means
        temp4$RESPONSE <- j
        temp4$YEN <- i
        temp4$DATA_TYPE <- "BLUE"
        blues_se <- rbind(blues_se, temp4)

        # Remove
        remove(temp2, temp3, temp4)
      } else if (unique(temp2$DATA_TYPE) == "BINOMIAL") {
        # Do nothing and just jam it into the blues dataframe
        temp3 <- data.frame(
          "rID" = temp2$rID,
          "predicted.value" = temp2$Y,
          "std.error" = NA,
          "RESPONSE" = j,
          "YEN" = i,
          "DATA_TYPE" = "BINOMIAL"
        )

        # Rbind
        blues_se <- rbind(blues_se, temp3)

        # Remove
        remove(temp2, temp3)
      }
    }
    # Remove
    remove(temp1)
  }

  # Write out
  write.csv(blues_se,
    "blues_single_environment.csv",
    row.names = FALSE
  )
} else {
  # Read in
  blues_se <- read.csv("blues_single_environment.csv")
}


## ----multi_env_adj--------------------------------------------------------------------------------------------------------
# If exists
if (!file.exists("blues_multiple_environment.csv")) {
  # Create object for ME blues
  blues_me <- blues_se |>
    dplyr::distinct(rID)

  # Save multienvironment models
  models_me <- list()

  # Loop
  for (i in unique(blues_se$RESPONSE)) {
    # Make message
    print_message <- paste("### Analyzing Trait =", i, "###")
    border <- paste(rep("#", nchar(print_message)), collapse = "")
    print(border)
    print(print_message)
    print(border)

    # Pull dataframe
    temp1 <- blues_se[blues_se$RESPONSE == i & blues_se$DATA_TYPE == "BLUE", ]

    # Make message
    print_message <- paste("### System Time =", Sys.time(), "###")
    border <- paste(rep("#", nchar(print_message)), collapse = "")
    print(border)
    print(print_message)
    print(border)

    # Fit mmes gBLUP for GWAS
    temp2 <- sommer::mmes(
      fixed = predicted.value ~ rID,
      random = ~YEN,
      rcov = ~units,
      W = diag(temp1$std.error),
      data = temp1
    )

    # Make message
    print_message <- paste("### System Time =", Sys.time(), "###")
    border <- paste(rep("#", nchar(print_message)), collapse = "")
    print(border)
    print(print_message)
    print(border)

    # Predict
    temp3 <- sommer::predict.mmes(temp2, D = "rID")$pvals

    # Rename
    colnames(temp3)[2:ncol(temp3)] <- paste(i, colnames(temp3)[2:ncol(temp3)], sep = "_")

    # Left join
    blues_me <- blues_me |> dplyr::left_join(temp3, by = "rID")

    # Save model
    models_me[[i]]$data <- temp1
    models_me[[i]]$model <- temp2

    # Remove
    remove(temp1, temp2, temp3)
  }

  # Write out adj means
  write.csv(blues_me,
    "blues_multiple_environment.csv",
    row.names = FALSE
  )

  # Make object for correlation visualization
  blues_me_corr <- blues_me |>
    dplyr::select(rID, GY_predicted.value, TW_predicted.value, Virus_predicted.value) |>
    dplyr::rename(
      `Grain Yield (Kg/Ha)` = GY_predicted.value,
      `Test Weight (Kg/Hl)` = TW_predicted.value,
      `Virus Ratings (Transformed)` = Virus_predicted.value
    ) |>
    tidyr::drop_na()

  # Now look at correlation among BLUEs
  jpeg(
    filename = "blues_multiple_environment_pairs_plotjpg",
    width = 8,
    height = 8,
    units = "in",
    res = 320
  )

  psych::pairs.panels(blues_me_corr[, 2:4],
    hist.col = "gray",
    stars = TRUE,
    lm = TRUE,
    line.col = "red",
    density = TRUE
  )

  dev.off()
} else {
  # Read data
  blues_me <- read.csv("blues_multiple_environment.csv")
}


## ----gwas_continuous------------------------------------------------------------------------------------------------------
# Read in data
geno <- gaston::read.vcf("virus_gwas_2025_production_final.vcf.gz", convert.chr = FALSE)

# Make a marker map
marker_map <- geno@snps |>
  dplyr::rename(
    chr_wheat = chr,
    snp_id = id
  ) |>
  dplyr::mutate(
    chr = as.numeric(as.factor(chr_wheat)),
    variant.id = c(1:nrow(geno@snps))
  ) |>
  dplyr::select(
    chr,
    chr_wheat,
    variant.id,
    snp_id,
    pos
  ) |>
  dplyr::mutate(
    chr = as.integer(chr),
    variant.id = as.integer(variant.id),
    pos = as.integer(pos)
  )

# Subset phenotypic infomration
pheno_file <- blues_me |>
  dplyr::filter(rID %in% geno@ped$id) |>
  dplyr::select(rID, GY_predicted.value, TW_predicted.value, Virus_predicted.value) |>
  dplyr::rename(
    GY = GY_predicted.value,
    TW = TW_predicted.value,
    Virus = Virus_predicted.value
  )

# Subset genotypic information
geno_file <- gaston::as.matrix(geno)
geno_file <- geno_file[rownames(geno_file) %in% pheno_file$rID, ]

# Now make a line map
line_map <- data.frame(
  rID = rownames(geno_file),
  scanID = 1:length(rownames(geno_file))
)

# Now bind the phenotype file with the scanID and genotype rownames
rownames(geno_file) <- line_map$scanID
pheno_file <- line_map |>
  dplyr::left_join(pheno_file, by = "rID")

# Make geontype data object for GENISIS
geno_genisis <- GWASTools::MatrixGenotypeReader(
  genotype = t(geno_file),
  snpID = marker_map$variant.id,
  chromosome = marker_map$chr,
  position = marker_map$pos,
  scanID = line_map$scanID
)

# Now make the genotype data
geno_data_genisis <- GWASTools::GenotypeData(geno_genisis)

# Calculate kinship matrix
kinship_file <- sommer::A.mat(geno_file - 1)

# Make message for system time
print_message <- paste("### System Time =", Sys.time(), "###")
border <- paste(rep("#", nchar(print_message)), collapse = "")
print(border)
print(print_message)
print(border)

# Do PCA
pcs_geno <- prcomp(geno_file,
  center = TRUE,
  scale. = TRUE
)

# Make message for system time
print_message <- paste("### System Time =", Sys.time(), "###")
border <- paste(rep("#", nchar(print_message)), collapse = "")
print(border)
print(print_message)
print(border)

# Pull percent variance
pve_geno <- summary(pcs_geno)$importance

# Create labels with variance explained
pve_labs <- paste0(
  "PC",
  1:ncol(pve_geno),
  " (",
  round(pve_geno[2, 1:ncol(pve_geno)] * 100, 2),
  "%)"
)

# Now pull PCs
pcs_geno <- pcs_geno$x

# Make dataframe
pcs_geno <- as.data.frame(pcs_geno)

# Generate and stitch plots
pcs_geno_image <- ggplot(pcs_geno, aes(PC1, PC2)) +
  geom_point() +
  theme_bw() +
  labs(x = pve_labs[1], y = pve_labs[2]) +
  ggplot(pcs_geno, aes(PC1, PC3)) +
  geom_point() +
  theme_bw() +
  labs(x = pve_labs[1], y = pve_labs[3]) +
  ggplot(pcs_geno, aes(PC2, PC3)) +
  geom_point() +
  theme_bw() +
  labs(x = pve_labs[2], y = pve_labs[3]) +
  patchwork::plot_annotation(title = "Principal Component Analysis of Genetic Data")

# Save image
ggsave(
  filename = "geno_pca_comparison_continuous.jpg",
  pcs_geno_image,
  width = 15,
  height = 5,
  units = "in",
  dpi = 320
)

# Make a new covar file for the final scanannot file
covar_file_scanAnnot <- pcs_geno |>
  tibble::rownames_to_column(var = "scanID")

# Pull first 10 PCs
covar_file_scanAnnot <- as.data.frame(covar_file_scanAnnot[, 1:11])
covar_file_scanAnnot$scanID <- as.integer(covar_file_scanAnnot$scanID)

# Make scanannot
scanAnnot <- pheno_file |>
  dplyr::left_join(covar_file_scanAnnot, by = "scanID") |>
  dplyr::select(scanID, PC1, PC2, PC3, GY, TW, Virus) |>
  GWASTools::ScanAnnotationDataFrame()

# Remove
remove(covar_file_scanAnnot)

# Plot function
plot_gwas <- function(df,
                      chr = "chr_wheat",
                      pos = "pos",
                      y = "neg_log_pval",
                      threshold = NULL,
                      title = NULL) {
  # We assume the input is sorted by chromosome and position
  plot_data <- data.frame(CHR = df[[chr]], POS = df[[pos]], Y = df[[y]])
  plot_data$CHR <- factor(plot_data$CHR, levels = unique(plot_data$CHR))

  # tapply gets the max POS per chromosome; cumsum creates the running total
  chr_max <- tapply(plot_data$POS, plot_data$CHR, max)
  offsets <- c(0, cumsum(as.numeric(chr_max))[-length(chr_max)])

  # Map offsets to the dataframe using the factor integer indices
  plot_data$BP_CUM <- plot_data$POS + offsets[as.numeric(plot_data$CHR)]

  # Calculate X-axis breaks (centers of each chromosome)
  axis_breaks <- tapply(plot_data$BP_CUM, plot_data$CHR, mean)

  # Generate plot
  plot_obj <- ggplot2::ggplot(plot_data, ggplot2::aes(x = BP_CUM, y = Y)) +
    ggplot2::geom_point(ggplot2::aes(color = CHR), alpha = 0.8, size = 1.5) +
    ggplot2::scale_color_manual(values = rep(c("grey", "black"), nlevels(plot_data$CHR))) +
    ggplot2::scale_x_continuous(label = names(axis_breaks), breaks = axis_breaks) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = "Chromosome", y = expression(-log[10](italic(P))))

  # if threshold != NULL
  if (!is.null(threshold)) {
    plot_obj <- plot_obj + ggplot2::geom_hline(
      yintercept = threshold,
      color = "red",
      linewidth = 1,
      linetype = "dashed"
    )
  }

  # if title != NULL
  if (!is.null(title)) {
    plot_obj <- plot_obj + ggplot2::labs(title = title)
  }

  # Return
  return(plot_obj)
}

# Make results object
GWAS_results_continuous <- c()

# Loop through traits
for (i in c("GY", "TW", "Virus")) {
  # Make message for system time
  print_message <- paste("### Analyzing Trait =", i, "###")
  border <- paste(rep("#", nchar(print_message)), collapse = "")
  print(border)
  print(print_message)
  print(border)

  # Make message for system time
  print_message <- paste("### System Time =", Sys.time(), "###")
  border <- paste(rep("#", nchar(print_message)), collapse = "")
  print(border)
  print(print_message)
  print(border)

  # Run Null Model for normal response
  temp1 <- GENESIS::fitNullModel(scanAnnot,
    outcome = i,
    covars = c("PC1", "PC2", "PC3"),
    cov.mat = kinship_file
  )

  # Genotype Block Iterator for GWAS
  temp2 <- GWASTools::GenotypeBlockIterator(geno_data_genisis)

  # Now run GWAS with null model
  temp3 <- suppressWarnings(GENESIS::assocTestSingle(temp2, null.model = temp1))

  # Now bind in true chromosome names and such
  temp4 <- marker_map |>
    dplyr::mutate(chr = as.character(chr)) |>
    dplyr::left_join(temp3, by = c("chr", "variant.id", "pos")) |>
    dplyr::mutate(
      neg_log_pval = -log10(Score.pval),
      trait = i
    ) |>
    dplyr::select(
      trait,
      chr,
      chr_wheat,
      variant.id,
      snp_id,
      pos,
      n.obs,
      freq,
      MAC,
      Score,
      Score.SE,
      Score.Stat,
      Score.pval,
      neg_log_pval,
      Est,
      Est.SE,
      PVE
    ) |>
    dplyr::mutate(break_bft = ifelse(neg_log_pval >= -log10(0.05 / nrow(temp3)),
      1,
      0
    ))

  # Now plot
  temp5 <- plot_gwas(temp4,
    threshold = -log10(0.05 / nrow(temp3))
  )

  # Now save image
  ggsave(
    filename = paste("gwas_manhattan_plot_", i, ".jpg", sep = ""),
    plot = temp5,
    width = 11,
    height = 8.5,
    units = "in",
    dpi = 320
  )

  # Now bind results
  GWAS_results_continuous <- rbind(GWAS_results_continuous, temp4)

  # Remove
  remove(temp1, temp2, temp3, temp4, temp5)

  # Make message for system time
  print_message <- paste("### System Time =", Sys.time(), "###")
  border <- paste(rep("#", nchar(print_message)), collapse = "")
  print(border)
  print(print_message)
  print(border)
}

# Write out results
write.csv(GWAS_results_continuous,
  "gwas_results_continuous.csv",
  row.names = FALSE
)


## ----gwas_binomial--------------------------------------------------------------------------------------------------------
# Pull observations for binomial ratings
binomial_ratings <- blues_se |>
  dplyr::filter(DATA_TYPE == "BINOMIAL") |>
  dplyr::select(YEN, rID, predicted.value) |>
  dplyr::rename(virus_binomial = predicted.value)

# Get the mode of each rating
binomial_ratings_consensus <- c()

# Loop through ratings
for (i in unique(binomial_ratings$rID)) {
  # Pull data
  temp1 <- binomial_ratings |>
    dplyr::filter(rID == i) |>
    dplyr::select(rID, virus_binomial)

  # Check if there is only one rating
  if (nrow(temp1) == 1) {
    # Bind
    binomial_ratings_consensus <- rbind(binomial_ratings_consensus, temp1)

    # Remove data
    remove(temp1)

    # Go to next in i
    next
  }

  # Now take table
  temp2 <- table(temp1$virus_binomial)

  # If table only has one rating
  if (length(temp2) == 1) {
    # Get single rating
    temp1 <- temp1 |>
      dplyr::distinct()

    # Bind
    binomial_ratings_consensus <- rbind(binomial_ratings_consensus, temp1)

    # Remove data
    remove(temp1, temp2)

    # Go to next in i
    next
  }

  # Take maximum
  temp3 <- names(temp2[temp2 == max(temp2)]) |> as.numeric()

  # If there is a tie
  if (length(temp3) > 1 | length(temp3) < 1) {
    # remove
    remove(temp1, temp2, temp3)

    # Skip
    next
  } else {
    # Make dataframe
    temp4 <- data.frame(
      rID = i,
      virus_binomial = temp3
    )

    # Bind
    binomial_ratings_consensus <- rbind(binomial_ratings_consensus, temp4)

    # Remove
    remove(temp1, temp2, temp3, temp4)
  }
}

# Make a marker map
marker_map <- geno@snps |>
  dplyr::rename(
    chr_wheat = chr,
    snp_id = id
  ) |>
  dplyr::mutate(
    chr = as.numeric(as.factor(chr_wheat)),
    variant.id = c(1:nrow(geno@snps))
  ) |>
  dplyr::select(
    chr,
    chr_wheat,
    variant.id,
    snp_id,
    pos
  ) |>
  dplyr::mutate(
    chr = as.integer(chr),
    variant.id = as.integer(variant.id),
    pos = as.integer(pos)
  )

# Subset phenotypic infomration
pheno_file <- binomial_ratings_consensus |>
  dplyr::filter(rID %in% geno@ped$id) |>
  dplyr::select(rID, virus_binomial) |>
  dplyr::rename(Virus = virus_binomial)

# Subset genotypic information
geno_file <- gaston::as.matrix(geno)
geno_file <- geno_file[rownames(geno_file) %in% pheno_file$rID, ]

# Now make a line map
line_map <- data.frame(
  rID = rownames(geno_file),
  scanID = 1:length(rownames(geno_file))
)

# Now bind the phenotype file with the scanID and genotype rownames
rownames(geno_file) <- line_map$scanID
pheno_file <- line_map |>
  dplyr::left_join(pheno_file, by = "rID")

# Make geontype data object for GENISIS
geno_genisis <- GWASTools::MatrixGenotypeReader(
  genotype = t(geno_file),
  snpID = marker_map$variant.id,
  chromosome = marker_map$chr,
  position = marker_map$pos,
  scanID = line_map$scanID
)

# Now make the genotype data
geno_data_genisis <- GWASTools::GenotypeData(geno_genisis)

# Calculate kinship matrix
kinship_file <- sommer::A.mat(geno_file - 1)

# Make message for system time
print_message <- paste("### System Time =", Sys.time(), "###")
border <- paste(rep("#", nchar(print_message)), collapse = "")
print(border)
print(print_message)
print(border)

# Do PCA
pcs_geno <- prcomp(geno_file,
  center = TRUE,
  scale. = TRUE
)

# Make message for system time
print_message <- paste("### System Time =", Sys.time(), "###")
border <- paste(rep("#", nchar(print_message)), collapse = "")
print(border)
print(print_message)
print(border)

# Pull percent variance
pve_geno <- summary(pcs_geno)$importance

# Create labels with variance explained
pve_labs <- paste0(
  "PC",
  1:ncol(pve_geno),
  " (",
  round(pve_geno[2, 1:ncol(pve_geno)] * 100, 2),
  "%)"
)

# Now pull PCs
pcs_geno <- pcs_geno$x

# Make dataframe
pcs_geno <- as.data.frame(pcs_geno)

# Generate and stitch plots
pcs_geno_image <- ggplot(pcs_geno, aes(PC1, PC2)) +
  geom_point() +
  theme_bw() +
  labs(x = pve_labs[1], y = pve_labs[2]) +
  ggplot(pcs_geno, aes(PC1, PC3)) +
  geom_point() +
  theme_bw() +
  labs(x = pve_labs[1], y = pve_labs[3]) +
  ggplot(pcs_geno, aes(PC2, PC3)) +
  geom_point() +
  theme_bw() +
  labs(x = pve_labs[2], y = pve_labs[3]) +
  patchwork::plot_annotation(title = "Principal Component Analysis of Genetic Data")

# Save image
ggsave(
  filename = "geno_pca_comparison_binomial.jpg",
  pcs_geno_image,
  width = 15,
  height = 5,
  units = "in",
  dpi = 320
)

# Make a new covar file for the final scanannot file
covar_file_scanAnnot <- pcs_geno |>
  tibble::rownames_to_column(var = "scanID")

# Pull first 10 PCs
covar_file_scanAnnot <- as.data.frame(covar_file_scanAnnot[, 1:11])
covar_file_scanAnnot$scanID <- as.integer(covar_file_scanAnnot$scanID)

# Make scanannot
scanAnnot <- pheno_file |>
  dplyr::left_join(covar_file_scanAnnot, by = "scanID") |>
  dplyr::select(scanID, PC1, PC2, PC3, Virus) |>
  GWASTools::ScanAnnotationDataFrame()

# Remove
remove(covar_file_scanAnnot)

# Make results object
GWAS_results_binomial <- c()

# Loop through traits
for (i in c("Virus")) {
  # Make message for system time
  print_message <- paste("### Analyzing Trait =", i, "###")
  border <- paste(rep("#", nchar(print_message)), collapse = "")
  print(border)
  print(print_message)
  print(border)

  # Make message for system time
  print_message <- paste("### System Time =", Sys.time(), "###")
  border <- paste(rep("#", nchar(print_message)), collapse = "")
  print(border)
  print(print_message)
  print(border)

  # Run Null Model for binomial response
  temp1 <- GENESIS::fitNullModel(scanAnnot,
    outcome = i,
    covars = c("PC1", "PC2", "PC3"),
    cov.mat = kinship_file,
    family = "binomial"
  )

  # Genotype Block Iterator for GWAS
  temp2 <- GWASTools::GenotypeBlockIterator(geno_data_genisis)

  # Now run GWAS with null model
  temp3 <- suppressWarnings(GENESIS::assocTestSingle(temp2, null.model = temp1))

  # Now bind in true chromosome names and such
  temp4 <- marker_map |>
    dplyr::mutate(chr = as.character(chr)) |>
    dplyr::left_join(temp3, by = c("chr", "variant.id", "pos")) |>
    dplyr::mutate(
      neg_log_pval = -log10(Score.pval),
      trait = i
    ) |>
    dplyr::select(
      trait,
      chr,
      chr_wheat,
      variant.id,
      snp_id,
      pos,
      n.obs,
      freq,
      MAC,
      Score,
      Score.SE,
      Score.Stat,
      Score.pval,
      neg_log_pval,
      Est,
      Est.SE,
      PVE
    ) |>
    dplyr::mutate(break_bft = ifelse(neg_log_pval >= -log10(0.05 / nrow(temp3)),
      1,
      0
    ))

  # Now plot
  temp5 <- plot_gwas(temp4,
    threshold = -log10(0.05 / nrow(temp3))
  )

  # Now save image
  ggsave(
    filename = paste("gwas_manhattan_plot_", i, "_binomial.jpg", sep = ""),
    plot = temp5,
    width = 11,
    height = 8.5,
    units = "in",
    dpi = 320
  )

  # Now bind results
  GWAS_results_binomial <- rbind(GWAS_results_binomial, temp4)

  # Remove
  remove(temp1, temp2, temp3, temp4, temp5)

  # Make message for system time
  print_message <- paste("### System Time =", Sys.time(), "###")
  border <- paste(rep("#", nchar(print_message)), collapse = "")
  print(border)
  print(print_message)
  print(border)
}

# Write out results
write.csv(GWAS_results_binomial,
  "gwas_results_binomial.csv",
  row.names = FALSE
)


## ----legacy_code, eval=FALSE, include=FALSE-------------------------------------------------------------------------------
# # # Read in data
# # geno <- gaston::read.vcf("virus_gwas_2025_production_final.vcf.gz", convert.chr = FALSE)
# #
# # #### statgenGWAS ####
# # # Make the marker map
# # marker_map <- geno@snps
# # rownames(marker_map) <- marker_map[, 2]
# # marker_map <- marker_map[, c(1, 4)]
# #
# # # Make the marker matrix
# # marker_mat <- gaston::as.matrix(geno)
# #
# # # Read in phenotypic file
# # pheno_file <- blues_me[, -grep("std.err", colnames(blues_me))]
# #
# # # Subset phenotyped individuals to only have individuals featured in marker matrix
# # pheno_file <- pheno_file[pheno_file$rID %in% rownames(marker_mat), ]
# #
# # # Subset marker matrix to only have phenotyped invdividuals
# # marker_mat <- marker_mat[rownames(marker_mat) %in% blues_me$rID, ]
# #
# # # Calculate kinship matrix
# # kinship_mat <- sommer::A.mat(marker_mat - 1)
# #
# # # Make message for system time
# # print_message <- paste("### System Time =", Sys.time(), "###")
# # border <- paste(rep("#", nchar(print_message)), collapse = "")
# # print(border)
# # print(print_message)
# # print(border)
# #
# # # Do PCA
# # pcs_geno <- prcomp(marker_mat,
# #   center = TRUE,
# #   scale. = TRUE
# # )
# #
# # # Make message for system time
# # print_message <- paste("### System Time =", Sys.time(), "###")
# # border <- paste(rep("#", nchar(print_message)), collapse = "")
# # print(border)
# # print(print_message)
# # print(border)
# #
# # # Pull percent variance
# # pve_geno <- summary(pcs_geno)$importance
# #
# # # Create labels with variance explained
# # pve_labs <- paste0(
# #   "PC",
# #   1:ncol(pve_geno),
# #   " (",
# #   round(pve_geno[2, 1:ncol(pve_geno)] * 100, 2),
# #   "%)"
# # )
# #
# # # Now pull PCs
# # pcs_geno <- pcs_geno$x
# #
# # # Make dataframe
# # pcs_geno <- as.data.frame(pcs_geno)
# #
# # # Generate and stitch plots
# # pcs_geno_image <- ggplot(pcs_geno, aes(PC1, PC2)) +
# #   geom_point() +
# #   theme_bw() +
# #   labs(x = pve_labs[1], y = pve_labs[2]) +
# #   ggplot(pcs_geno, aes(PC1, PC3)) +
# #   geom_point() +
# #   theme_bw() +
# #   labs(x = pve_labs[1], y = pve_labs[3]) +
# #   ggplot(pcs_geno, aes(PC2, PC3)) +
# #   geom_point() +
# #   theme_bw() +
# #   labs(x = pve_labs[2], y = pve_labs[3]) +
# #   patchwork::plot_annotation(title = "Principal Component Analysis of Genetic Data")
# #
# # # Save image
# # ggsave(
# #   filename = "geno_pca_comparison.jpg",
# #   pcs_geno_image,
# #   width = 15,
# #   height = 5,
# #   units = "in",
# #   dpi = 320
# # )
# #
# # # Make a covariate file with just lines in phenotype (just to make sure)
# # covar_file <- pcs_geno[rownames(pcs_geno) %in% pheno_file$rID, ]
# #
# # # Make sure first column is titled "genotype"
# # colnames(pheno_file)[1] <- "genotype"
# #
# # # Make gdata object for statgenGWAS
# # gdata_obj <- statgenGWAS::createGData(
# #   geno = marker_mat,
# #   map = marker_map,
# #   kin = kinship_mat,
# #   pheno = pheno_file,
# #   covar = covar_file
# # )
# #
# # # Make message for system time
# # print_message <- paste("### System Time =", Sys.time(), "###")
# # border <- paste(rep("#", nchar(print_message)), collapse = "")
# # print(border)
# # print(print_message)
# # print(border)
# #
# # # Run GWAS on all traits
# # GWAS_Results <- statgenGWAS::runSingleTraitGwas(
# #   gData = gdata_obj,
# #   traits = colnames(pheno_file)[2:ncol(pheno_file)],
# #   covar = c("PC1", "PC2", "PC3"),
# #   kin = gdata_obj$kinship
# # )
# #
# # # Make message for system time
# # print_message <- paste("### System Time =", Sys.time(), "###")
# # border <- paste(rep("#", nchar(print_message)), collapse = "")
# # print(border)
# # print(print_message)
# # print(border)
# #
# # # Quick plot
# # gy_gwas_img <- plot(GWAS_Results,
# #   plotType = "manhattan",
# #   trait = "GY_predicted.value",
# #   title = "Grain Yield (Kg/Ha)"
# # )
# # tw_gwas_img <- plot(GWAS_Results,
# #   plotType = "manhattan",
# #   trait = "TW_predicted.value",
# #   title = "Test Weight (Kg/Hl)"
# # )
# # virus_gwas_img <- plot(GWAS_Results,
# #   plotType = "manhattan",
# #   trait = "Virus_predicted.value",
# #   title = "Virus (Transformed)"
# # )
# #
# # # Write images
# # ggsave("gwas_image_gy.jpg",
# #   plot = gy_gwas_img,
# #   width = 20,
# #   height = 10,
# #   units = "in",
# #   dpi = 320
# # )
# # ggsave("gwas_image_tw.jpg",
# #   plot = tw_gwas_img,
# #   width = 20,
# #   height = 10,
# #   units = "in",
# #   dpi = 320
# # )
# # ggsave("gwas_image_virus.jpg",
# #   plot = virus_gwas_img,
# #   width = 20,
# #   height = 10,
# #   units = "in",
# #   dpi = 320
# # )
# #
# # # Write out results
# # write.csv(GWAS_Results$GWAResult$pheno_file,
# #   "gwas_results_full.csv",
# #   row.names = FALSE
# # )
#
# # #### rrBLUP ####
# # # Make marker matrix
# # marker_map <- geno@snps
# # marker_map <- marker_map[, c(2,1,4)]
# # marker_mat <- t(gaston::as.matrix(geno)-1) |>
# #   as.data.frame() |>
# #   tibble::rownames_to_column(var = "id") |>
# #   as.data.frame()
# # marker_mat <- marker_map |>
# #   dplyr::left_join(marker_mat, by = "id")
# #
# # # Make phenotype file
# # pheno_file <- blues_me[, -grep("std.err", colnames(blues_me))]
# # pheno_file <- pheno_file[pheno_file$rID %in% colnames(marker_mat), ]
# #
# # # Make sure that the individuals who are in marker matrix are in phenotypic matrix
# # lines <- colnames(marker_mat)[colnames(marker_mat) %in% pheno_file$rID]
# # marker_mat <- marker_mat[,c(colnames(marker_mat)[1:3],
# #                             colnames(marker_mat)[colnames(marker_mat) %in% lines])]
# #
# # # Plot function
# # plot_gwas <- function(df, chr = "chr", pos = "position", y = "-log(p)") {
# #
# #   # Create a local subset with standardized names to avoid tidy-eval dependencies
# #   # We assume the input is sorted by chromosome and position
# #   plot_data <- data.frame(CHR = df[[chr]], POS = df[[pos]], Y = df[[y]])
# #   plot_data$CHR <- factor(plot_data$CHR, levels = unique(plot_data$CHR))
# #
# #   # Calculate cumulative position offsets using base R
# #   # tapply gets the max POS per chromosome; cumsum creates the running total
# #   chr_max <- tapply(plot_data$POS, plot_data$CHR, max)
# #   offsets <- c(0, cumsum(as.numeric(chr_max))[-length(chr_max)])
# #
# #   # Map offsets to the dataframe using the factor integer indices
# #   plot_data$BP_CUM <- plot_data$POS + offsets[as.numeric(plot_data$CHR)]
# #
# #   # Calculate X-axis breaks (centers of each chromosome)
# #   axis_breaks <- tapply(plot_data$BP_CUM, plot_data$CHR, mean)
# #
# #   # Generate plot
# #   ggplot2::ggplot(plot_data, ggplot2::aes(x = BP_CUM, y = Y)) +
# #     ggplot2::geom_point(ggplot2::aes(color = CHR), alpha = 0.8, size = 1.5) +
# #     ggplot2::scale_color_manual(values = rep(c("grey30", "navy"), nlevels(plot_data$CHR))) +
# #     ggplot2::scale_x_continuous(label = names(axis_breaks), breaks = axis_breaks) +
# #     ggplot2::theme_minimal() +
# #     ggplot2::theme(legend.position = "none",
# #                   panel.grid.major.x = ggplot2::element_blank(),
# #                   panel.grid.minor.x = ggplot2::element_blank()) +
# #     ggplot2::labs(x = "Chromosome", y = expression(-log[10](italic(P))))
# # }
# #
# # # Run GWAS
# # GWAS_rrBLUP_Full <- rrBLUP::GWAS(pheno = pheno_file,
# #                                  geno = marker_mat,
# #                                  n.PC = 3,
# #                                  plot = FALSE)
# #
# # # Plot
# # plot_gwas(GWAS_rrBLUP_Full, chr = "chr", pos = "pos", y = "GY_predicted.value")
# #
# # # Run GWAS
# # GWAS_rrBLUP_Sub <- rrBLUP::GWAS(pheno = pheno_file,
# #                                 geno = marker_mat[marker_mat$chr %in% c("4A"),],
# #                                 n.PC = 3,
# #                                 plot = FALSE)
# # #### Save ####
# # save.image(file = "Final_Image.RData")
