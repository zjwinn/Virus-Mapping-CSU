## ----setup, include=FALSE----------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----single_env_adj----------------------------------------------------------------------------------------------
# Ensure BiocManager is installed, as GWASTools and GENESIS are Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Define the required packages
target_packages <- c("GWASTools", "GENESIS")

# Identify which packages are not currently installed
missing_packages <- target_packages[!(target_packages %in% installed.packages()[, "Package"])]

# Install the missing packages if any
if (length(missing_packages) > 0) {
  cat("Installing missing packages: ", paste(missing_packages, collapse = ", "), "\n")
  # Setting update = FALSE prevents it from prompting to update all other Bioc packages,
  # ask = FALSE prevents interactive prompts
  BiocManager::install(missing_packages, update = FALSE, ask = FALSE)
  cat("Installation complete.\n")
} else {
  cat("Success: GWASTools and GENESIS are already installed.\n")
}

# Make a directory
directory <- file.path(
  getwd(),
  "single_environment_adjustment"
)
dir.create(directory, showWarnings = FALSE)

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

# Pre yellowing disease (wsmv)
pre_yd <- pheno |>
  dplyr::group_by(YEN) |>
  dplyr::summarize(
    Traits = paste(unique(RESPONSE), collapse = ", "),
    .groups = "drop"
  ) |>
  dplyr::filter(Traits %in% Traits[grep("WSMV", Traits)]) |>
  tidyr::separate(YEN, sep = ":", into = c("Year", "Env", "Nusery"), remove = FALSE) |>
  dplyr::mutate(Year = as.numeric(Year)) |>
  dplyr::filter(Year <= 2019) |>
  dplyr::select(YEN)
pre_yd <- pre_yd$YEN

# Post yellowing disease (trimv)
post_yd <- pheno |>
  dplyr::group_by(YEN) |>
  dplyr::summarize(
    Traits = paste(unique(RESPONSE), collapse = ", "),
    .groups = "drop"
  ) |>
  dplyr::filter(Traits %in% Traits[grep("Virus", Traits)]) |>
  tidyr::separate(YEN, sep = ":", into = c("Year", "Env", "Nusery"), remove = FALSE) |>
  dplyr::mutate(Year = as.numeric(Year)) |>
  dplyr::filter(Year >= 2021) |>
  dplyr::select(YEN)
post_yd <- post_yd$YEN

# Get the total number of environments we want
yen_total <- unique(c(pre_yd, post_yd))

# Filter data
pheno <- pheno |>
  dplyr::filter(YEN %in% yen_total)

# Get summary
pheno_summary <- pheno |>
  dplyr::group_by(YEN, RESPONSE) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
  tidyr::pivot_wider(names_from = RESPONSE, values_from = n, values_fill = 0)

# Pre summary
pheno_summary_pre <- pheno_summary |>
  dplyr::filter(YEN %in% pre_yd) |>
  dplyr::select(YEN, GY, TW, WSMV) |>
  tidyr::separate(YEN, into = c("Year", "Location", "Nursery"), sep = ":") |>
  dplyr::mutate(Location = paste(Location, ", CO", sep = "")) |>
  dplyr::rename(
    `Grain Yield` = GY,
    `Test Weight` = TW,
    `WSMV Ratings` = WSMV
  )

# Post summaries
pheno_summary_post <- pheno_summary |>
  dplyr::filter(YEN %in% post_yd) |>
  dplyr::select(YEN, GY, TW, Virus) |>
  tidyr::separate(YEN, into = c("Year", "Location", "Nursery"), sep = ":") |>
  dplyr::mutate(Location = paste(Location, ", CO", sep = "")) |>
  dplyr::rename(
    `Grain Yield` = GY,
    `Test Weight` = TW,
    `Novel Disease Ratings` = Virus
  )

# Get summary of number of obs
colSums(pheno_summary_pre[, 4:ncol(pheno_summary_pre)], na.rm = TRUE)
colSums(pheno_summary_post[, 4:ncol(pheno_summary_post)], na.rm = TRUE)

# Get count of ordinal observations
pheno_summary_post_binomial <- pheno_summary |>
  dplyr::filter(YEN %in% YEN[grep("2025:Julesburg", YEN)]) |>
  dplyr::select(YEN, Virus) |>
  tidyr::separate(YEN, into = c("Year", "Location", "Nursery"), sep = ":") |>
  dplyr::mutate(Location = paste(Location, ", CO", sep = "")) |>
  dplyr::rename(`Virus Ratings` = Virus)

# Get some numbers for paper
sum(pheno_summary_post_binomial$`Virus Ratings`)
colSums(pheno_summary_post[, 4:ncol(pheno_summary_post)], na.rm = TRUE)[3] - sum(pheno_summary_post_binomial$`Virus Ratings`)

# Now get table 1
table_1 <- pheno_summary |>
  tidyr::separate(YEN, into = c("Year", "Location", "Nursery"), sep = ":", remove = FALSE) |>
  dplyr::group_by(Year, Location) |>
  dplyr::summarise(
    `Grain Yield` = sum(GY, na.rm = TRUE),
    `Test Weight` = sum(TW, na.rm = TRUE),
    `WSMV Visual Ratings - Ordinal` = sum(WSMV, na.rm = TRUE),
    `NDVR - Ordinal` = sum(Virus[!grepl("2025:Julesburg", YEN)], na.rm = TRUE),
    `NDVR - Binomial` = sum(Virus[grepl("2025:Julesburg", YEN)], na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(Year, Location)

# Write out table 1
write.csv(table_1,
  file.path(directory, "Table_1.csv"),
  row.names = FALSE
)

# Check if the analysis needs to be ran
if (!file.exists(file.path(directory, "blues_single_environment.csv"))) {
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
        temp3 <- try(
          sommer::mmes(
            fixed = Y ~ rID,
            random = ~ sommer::spl2Dc(ROW, COLUMN),
            rcov = ~units,
            data = temp2,
            dateWarning = FALSE,
            henderson = FALSE
          )
        )

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
        temp3 <- try(
          sommer::mmes(
            fixed = Y_t ~ rID,
            random = ~ sommer::spl2Dc(ROW, COLUMN),
            rcov = ~units,
            data = temp2,
            dateWarning = FALSE,
            henderson = FALSE
          )
        )

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
    file.path(directory, "blues_single_environment.csv"),
    row.names = FALSE
  )
} else {
  # Read in
  blues_se <- read.csv(file.path(directory, "blues_single_environment.csv"))
}


## ----multi_env_adj-----------------------------------------------------------------------------------------------
# Make a directory
directory <- file.path(
  getwd(),
  "multiple_environment_adjustment"
)
dir.create(directory, showWarnings = FALSE)

# If exists
if (!file.exists(file.path(directory, "blues_multiple_environment_pre_yd.csv")) |
  !file.exists(file.path(directory, "blues_multiple_environment_post_yd.csv"))) {
  #### Pre Yellowing Disease ####

  # Create object for ME blues
  blues_me_pre <- blues_se |>
    dplyr::filter(YEN %in% pre_yd) |>
    dplyr::distinct(rID)

  # Save multienvironment models
  models_me_pre <- list()

  # Divide out pre yd
  blues_se_pre_yd <- blues_se |>
    dplyr::filter(YEN %in% pre_yd)

  # Loop
  for (i in unique(blues_se_pre_yd$RESPONSE)) {
    # Make message
    print_message <- paste("### Analyzing Trait =", i, "###")
    border <- paste(rep("#", nchar(print_message)), collapse = "")
    print(border)
    print(print_message)
    print(border)

    # Pull dataframe
    temp1 <- blues_se_pre_yd[blues_se_pre_yd$RESPONSE == i & blues_se_pre_yd$DATA_TYPE == "BLUE", ]

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
      W = diag(1 / temp1$std.error^2), ,
      data = temp1,
      dateWarning = FALSE,
      henderson = FALSE
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
    blues_me_pre <- blues_me_pre |> dplyr::left_join(temp3, by = "rID")

    # Save model
    models_me_pre[[i]]$data <- temp1
    models_me_pre[[i]]$model <- temp2

    # Remove
    remove(temp1, temp2, temp3)
  }

  # Write out adj means
  write.csv(
    blues_me_pre,
    file.path(directory, "blues_multiple_environment_pre_yd.csv"),
    row.names = FALSE
  )

  # Make object for correlation visualization
  blues_me_corr <- blues_me_pre |>
    dplyr::select(
      rID,
      GY_predicted.value,
      TW_predicted.value,
      WSMV_predicted.value
    ) |>
    dplyr::rename(
      `Grain Yield (Kg/Ha)` = GY_predicted.value,
      `Test Weight (Kg/Hl)` = TW_predicted.value,
      `WSMV (Transformed)` = WSMV_predicted.value
      # `Virus Ratings (Transformed)` = Virus_predicted.value
    ) |>
    tidyr::drop_na()

  # Now look at correlation among BLUEs
  jpeg(
    filename = file.path(directory, "blues_multiple_environment_pre_yd_pairs_plot.jpg"),
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

  #### Post Yellowing Disease ####

  # Create object for ME blues
  blues_me_post <- blues_se |>
    dplyr::filter(YEN %in% post_yd) |>
    dplyr::distinct(rID)

  # Save multienvironment models
  models_me_post <- list()

  # Divide out pre yd
  blues_se_post_yd <- blues_se |>
    dplyr::filter(YEN %in% post_yd)

  # Loop
  for (i in unique(blues_se_post_yd$RESPONSE)) {
    # Make message
    print_message <- paste("### Analyzing Trait =", i, "###")
    border <- paste(rep("#", nchar(print_message)), collapse = "")
    print(border)
    print(print_message)
    print(border)

    # Pull dataframe
    temp1 <- blues_se_post_yd[blues_se_post_yd$RESPONSE == i & blues_se_post_yd$DATA_TYPE == "BLUE", ]

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
      W = diag(1 / temp1$std.error^2),
      data = temp1,
      dateWarning = FALSE,
      henderson = FALSE
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
    blues_me_post <- blues_me_post |> dplyr::left_join(temp3, by = "rID")

    # Save model
    models_me_post[[i]]$data <- temp1
    models_me_post[[i]]$model <- temp2

    # Remove
    remove(temp1, temp2, temp3)
  }

  # Write out adj means
  write.csv(blues_me_post,
    file.path(directory, "blues_multiple_environment_post_yd.csv"),
    row.names = FALSE
  )

  # Make object for correlation visualization
  blues_me_corr <- blues_me_post |>
    dplyr::select(rID, GY_predicted.value, TW_predicted.value, Virus_predicted.value) |>
    dplyr::rename(
      `Grain Yield (Kg/Ha)` = GY_predicted.value,
      `Test Weight (Kg/Hl)` = TW_predicted.value,
      # `WSMV (Transformed)` = WSMV_predicted.value
      `Virus Ratings (Transformed)` = Virus_predicted.value
    ) |>
    tidyr::drop_na()

  # Now look at correlation among BLUEs
  jpeg(
    filename = file.path(directory, "blues_multiple_environment_post_yd_pairs_plot.jpg"),
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

  # Calculate Table 2
  table_2 <- blues_me_pre |>
    dplyr::select(2, 4, 6) |>
    dplyr::rename(
      "Grain Yield (Kg/Ha)" = 1,
      "Test Weight (Kg/Hl)" = 2,
      "WSMV Visual Ratings (Transformed)" = 3
    ) |>
    tidyr::pivot_longer(c("Grain Yield (Kg/Ha)", "Test Weight (Kg/Hl)", "WSMV Visual Ratings (Transformed)")) |>
    dplyr::rename("Trait" = name) |>
    dplyr::group_by(Trait) |>
    dplyr::summarise(
      Min = min(value, na.rm = TRUE),
      Mean = mean(value, na.rm = TRUE),
      Max = max(value, na.rm = TRUE),
      SD = sd(value, na.rm = TRUE)
    ) |>
    dplyr::mutate(`Data Set` = "PreDCS") |>
    dplyr::select(`Data Set`, dplyr::everything())

  # Temp for table 2
  temp <- blues_me_post |>
    dplyr::select(2, 4, 6) |>
    dplyr::rename(
      "Grain Yield (Kg/Ha)" = 1,
      "Test Weight (Kg/Hl)" = 2,
      "Viral Visual Ratings (Transformed)" = 3
    ) |>
    tidyr::pivot_longer(c("Grain Yield (Kg/Ha)", "Test Weight (Kg/Hl)", "Viral Visual Ratings (Transformed)")) |>
    dplyr::rename("Trait" = name) |>
    dplyr::group_by(Trait) |>
    dplyr::summarise(
      Min = min(value, na.rm = TRUE),
      Mean = mean(value, na.rm = TRUE),
      Max = max(value, na.rm = TRUE),
      SD = sd(value, na.rm = TRUE)
    ) |>
    dplyr::mutate(`Data Set` = "PostDCS") |>
    dplyr::select(`Data Set`, dplyr::everything())

  # Rbind
  table_2 <- rbind(table_2, temp)

  # Order
  table_2 <- table_2[c(1, 4, 2, 5, 3, 6), ]

  # Write out
  write.csv(
    table_2,
    file.path(directory, "Table_2.csv"),
    row.names = FALSE
  )

  # Save image
  save.image(file.path(directory, "multi-environmental-models.RData"))
} else {
  # Read data
  blues_me_pre <- read.csv(file.path(directory, "blues_multiple_environment_pre_yd.csv"))
  blues_me_post <- read.csv(file.path(directory, "blues_multiple_environment_post_yd.csv"))
}


## ----gwas_continuous---------------------------------------------------------------------------------------------
# Make a directory
directory <- file.path(
  getwd(),
  "PreDCS_GWAS"
)
dir.create(directory, showWarnings = FALSE)

#### GWAS Plotting Function ####

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

#### Pre Yellowing Disease ####

# Read in data
geno <- gaston::read.vcf("virus_mapping_2026_production_final.vcf.gz", convert.chr = FALSE)

# Reassign chromosome names
geno@snps$chr <- gsub("Chr", "", geno@snps$chr)

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
pheno_file <- blues_me_pre |>
  dplyr::filter(rID %in% geno@ped$id) |>
  dplyr::select(rID, GY_predicted.value, TW_predicted.value, WSMV_predicted.value) |>
  dplyr::rename(
    GY = GY_predicted.value,
    TW = TW_predicted.value,
    WSMV = WSMV_predicted.value
  )

# Get numbers
blues_me_pre_summary <- pheno_file |>
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(!is.na(.)))) |>
  dplyr::select(-rID)

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
pcs_geno <- prcomp(geno_file)

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
pcs_geno_image <- ggplot2::ggplot(pcs_geno, ggplot2::aes(PC1, PC2)) +
  ggplot2::geom_point() +
  ggplot2::theme_bw() +
  ggplot2::labs(x = pve_labs[1], y = pve_labs[2]) +
  ggplot2::ggplot(pcs_geno, ggplot2::aes(PC1, PC3)) +
  ggplot2::geom_point() +
  ggplot2::theme_bw() +
  ggplot2::labs(x = pve_labs[1], y = pve_labs[3]) +
  ggplot2::ggplot(pcs_geno, ggplot2::aes(PC2, PC3)) +
  ggplot2::geom_point() +
  ggplot2::theme_bw() +
  ggplot2::labs(x = pve_labs[2], y = pve_labs[3]) +
  patchwork::plot_annotation(title = "Principal Component Analysis of Genetic Data")

# Save image
ggplot2::ggsave(
  filename = file.path(directory, "geno_pca_comparison_continuous_pre_dcs.jpg"),
  pcs_geno_image,
  width = 15,
  height = 5,
  units = "in",
  dpi = 320
)

# Make a new covar file for the final scanannot file
covar_file_scanAnnot <- pcs_geno |>
  tibble::rownames_to_column(var = "scanID")

# Make RDS of information for image
for_rds_obj <- list(
  pcs_geno = pcs_geno,
  pve_geno = pve_geno,
  line_map = line_map
)

# Save RDS
saveRDS(
  for_rds_obj,
  file = file.path(
    directory,
    "pca_results.RDS"
  )
)

# Remove RDS
remove(for_rds_obj)

# Pull first 10 PCs
covar_file_scanAnnot <- as.data.frame(covar_file_scanAnnot[, 1:11])
covar_file_scanAnnot$scanID <- as.integer(covar_file_scanAnnot$scanID)

# Make scanannot
scanAnnot <- pheno_file |>
  dplyr::left_join(covar_file_scanAnnot, by = "scanID") |>
  dplyr::select(scanID, PC1, PC2, PC3, GY, TW, WSMV) |>
  GWASTools::ScanAnnotationDataFrame()

# Remove
remove(covar_file_scanAnnot)

# Make results object
GWAS_results_continuous_pre_yd <- c()

# Loop through traits
for (i in c("GY", "TW", "WSMV")) {
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
  temp2 <- GWASTools::GenotypeBlockIterator(geno_data_genisis, snpBlock = 1000)

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
    dplyr::mutate(
      break_bft = ifelse(neg_log_pval >= -log10(0.05 / nrow(temp3)),
        1,
        0
      )
    )

  # Now plot
  temp5 <- plot_gwas(
    temp4,
    threshold = -log10(0.05 / nrow(temp3))
  )

  # Now save image
  ggplot2::ggsave(
    filename = file.path(directory, paste("gwas_manhattan_plot_pre_yd_", i, ".jpg", sep = "")),
    plot = temp5,
    width = 11,
    height = 8.5,
    units = "in",
    dpi = 320
  )

  # Now bind results
  GWAS_results_continuous_pre_yd <- rbind(GWAS_results_continuous_pre_yd, temp4)

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
write.csv(
  GWAS_results_continuous_pre_yd,
  file.path(directory, "gwas_results_continuous_pre_yd.csv"),
  row.names = FALSE
)

#### Post Yellowing Disease ####
# Make a directory
directory <- file.path(
  getwd(),
  "PostDCS_GWAS"
)
dir.create(directory, showWarnings = FALSE)

# Subset phenotypic infomration
pheno_file <- blues_me_post |>
  dplyr::filter(rID %in% geno@ped$id) |>
  dplyr::select(rID, GY_predicted.value, TW_predicted.value, Virus_predicted.value) |>
  dplyr::rename(
    GY = GY_predicted.value,
    TW = TW_predicted.value,
    Virus = Virus_predicted.value
  )

# Get numbers
blues_me_post_summary <- pheno_file |>
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(!is.na(.)))) |>
  dplyr::select(-rID)

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
pcs_geno <- prcomp(geno_file)

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
pcs_geno_image <- ggplot2::ggplot(pcs_geno, ggplot2::aes(PC1, PC2)) +
  ggplot2::geom_point() +
  ggplot2::theme_bw() +
  ggplot2::labs(x = pve_labs[1], y = pve_labs[2]) +
  ggplot2::ggplot(pcs_geno, ggplot2::aes(PC1, PC3)) +
  ggplot2::geom_point() +
  ggplot2::theme_bw() +
  ggplot2::labs(x = pve_labs[1], y = pve_labs[3]) +
  ggplot2::ggplot(pcs_geno, ggplot2::aes(PC2, PC3)) +
  ggplot2::geom_point() +
  ggplot2::theme_bw() +
  ggplot2::labs(x = pve_labs[2], y = pve_labs[3]) +
  patchwork::plot_annotation(title = "Principal Component Analysis of Genetic Data")

# Save image
ggplot2::ggsave(
  filename = file.path(directory, "geno_pca_comparison_continuous_post_yd.jpg"),
  pcs_geno_image,
  width = 15,
  height = 5,
  units = "in",
  dpi = 320
)

# Make RDS of information for image
for_rds_obj <- list(
  pcs_geno = pcs_geno,
  pve_geno = pve_geno,
  line_map = line_map
)

# Save RDS
saveRDS(
  for_rds_obj,
  file = file.path(
    directory,
    "pca_results.RDS"
  )
)

# Remove RDS
remove(for_rds_obj)

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

# Make results object
GWAS_results_continuous_post_yd <- c()

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
  temp2 <- GWASTools::GenotypeBlockIterator(geno_data_genisis, snpBlock = 1000)

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
  temp5 <- plot_gwas(
    temp4,
    threshold = -log10(0.05 / nrow(temp3))
  )

  # Now save image
  ggplot2::ggsave(
    filename = file.path(directory, paste("gwas_manhattan_plot_post_yd_", i, ".jpg", sep = "")),
    plot = temp5,
    width = 11,
    height = 8.5,
    units = "in",
    dpi = 320
  )

  # Now bind results
  GWAS_results_continuous_post_yd <- rbind(GWAS_results_continuous_post_yd, temp4)

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
write.csv(GWAS_results_continuous_post_yd,
  file.path(directory, "gwas_results_continuous_post_yd.csv"),
  row.names = FALSE
)


## ----gwas_binomial-----------------------------------------------------------------------------------------------
# Make a directory
directory <- file.path(
  getwd(),
  "PostDCS_GWAS_binomial"
)
dir.create(directory, showWarnings = FALSE)

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

binomial_ratings_consensus_summary <- binomial_ratings_consensus |>
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(!is.na(.)))) |>
  dplyr::select(-rID)

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
pcs_geno <- prcomp(geno_file)

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
pcs_geno_image <- ggplot2::ggplot(pcs_geno, ggplot2::aes(PC1, PC2)) +
  ggplot2::geom_point() +
  ggplot2::theme_bw() +
  ggplot2::labs(x = pve_labs[1], y = pve_labs[2]) +
  ggplot2::ggplot(pcs_geno, ggplot2::aes(PC1, PC3)) +
  ggplot2::geom_point() +
  ggplot2::theme_bw() +
  ggplot2::labs(x = pve_labs[1], y = pve_labs[3]) +
  ggplot2::ggplot(pcs_geno, ggplot2::aes(PC2, PC3)) +
  ggplot2::geom_point() +
  ggplot2::theme_bw() +
  ggplot2::labs(x = pve_labs[2], y = pve_labs[3]) +
  patchwork::plot_annotation(title = "Principal Component Analysis of Genetic Data")

# Save image
ggplot2::ggsave(
  filename = file.path(directory, "geno_pca_comparison_binomial.jpg"),
  pcs_geno_image,
  width = 15,
  height = 5,
  units = "in",
  dpi = 320
)

# Make RDS of information for image
for_rds_obj <- list(
  pcs_geno = pcs_geno,
  pve_geno = pve_geno,
  line_map = line_map
)

# Save RDS
saveRDS(
  for_rds_obj,
  file = file.path(
    directory,
    "pca_results.RDS"
  )
)

# Remove RDS
remove(for_rds_obj)

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
  temp2 <- GWASTools::GenotypeBlockIterator(geno_data_genisis, snpBlock = 1000)

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
  ggplot2::ggsave(
    filename = file.path(directory, paste("gwas_manhattan_plot_", i, "_binomial.jpg", sep = "")),
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
write.csv(
  GWAS_results_binomial,
  file.path(directory, "gwas_results_binomial.csv"),
  row.names = FALSE
)

