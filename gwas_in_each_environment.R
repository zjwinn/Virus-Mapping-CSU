# Plot function
plot_gwas <- function(
  df,
  chr = "chr_wheat",
  pos = "pos",
  y = "neg_log_pval",
  threshold = NULL,
  title = NULL
) {
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

# Read in BLUEs of single environments
blues_se <- read.csv(
  file.path(
    getwd(),
    "single_environment_adjustment",
    "blues_single_environment.csv"
  )
)

# Read in genotype file
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

# get the other columns we need
blues_se <- blues_se |>
  tidyr::separate(col = "YEN", into = c("YEAR", "EXPT", "NURNAME"), sep = ":", remove = FALSE) |>
  dplyr::mutate(YE = paste(YEAR, EXPT, sep = ":"))


# Define thresholds
maf_threshold <- 0.05
hz_threshold <- 0.10

# Make list_dir
list_path <- file.path(
  getwd(),
  "single_environment_GWAS"
)

# Make dir
dir.create(
  list_path,
  showWarnings = FALSE
)

# Now permutation per YE
for (i in unique(blues_se$YE)) {
  # Make message for system time
  print_message <- paste("### YEAR:EXPT =", i, "###")
  border <- paste(rep("#", nchar(print_message)), collapse = "")
  print(border)
  print(print_message)
  print(border)

  # Pull data
  temp_pheno <- blues_se |>
    dplyr::filter(YE == i) |>
    dplyr::filter(DATA_TYPE == "BLUE") |>
    tidyr::drop_na(RESPONSE)

  # Make temp file path
  temp_filepath <- file.path(list_path, gsub(":", "_", i))

  # Make a directory
  dir.create(
    temp_filepath,
    showWarnings = FALSE
  )

  # Now perm per trait
  for (j in unique(temp_pheno$RESPONSE)) {
    # Make message for system time
    print_message <- paste("### Response =", j, "###")
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

    # Pull data
    temp_pheno_sub <- temp_pheno |>
      dplyr::filter(RESPONSE == j)

    # If there is only one nursery in the environment
    if (length(unique(temp_pheno_sub$YEN)) == 1) {
      print(paste0("There is only one nusery for ", j, " in ", i, ". Skipping multi-nursery adjustment..."))

      # Make an adjusted means object
      temp_adj_means <- temp_pheno_sub |>
        dplyr::select(rID, predicted.value, std.error)
    } else {
      # Now run model
      temp_fit <- try(
        sommer::mmes(
          fixed = predicted.value ~ rID,
          random = ~YE,
          rcov = ~units,
          W = diag(1 / (temp_pheno_sub$std.error^2)),
          data = temp_pheno_sub,
          dateWarning = FALSE,
          henderson = FALSE
        )
      )

      # If try fail
      if (class(temp_fit) == "try-error") {
        # Say something
        print(paste("Multi-environment MLM of", j, "in", i, "failed. Skipping..."))

        # Remove
        remove(
          temp_pheno_sub,
          temp_fit
        )

        # Next
        next()
      }

      # Pull adjusted means
      temp_adj_means <- sommer::predict.mmes(temp_fit, D = "rID")$pvals

      # Remove the temp_fit object
      remove(temp_fit)
    }

    # Pull geno
    temp_geno <- gaston::as.matrix(geno)
    temp_geno <- temp_geno[rownames(temp_geno) %in% temp_adj_means$rID, ]

    # Calculate alternate allele frequency
    alt_freq <- colMeans(temp_geno, na.rm = TRUE) / 2

    # Calculate minor allele frequency
    maf <- pmin(alt_freq, 1 - alt_freq)

    # Calculate observed heterozygosity
    hz <- colMeans(temp_geno == 1, na.rm = TRUE)

    # Identify markers meeting the criteria
    keep_markers <- (maf >= maf_threshold) & (hz <= hz_threshold)

    # Subset the genotype matrix
    temp_geno <- temp_geno[, keep_markers]

    # Remove
    remove(
      alt_freq,
      maf,
      hz,
      keep_markers
    )

    # Make temp line map
    temp_line_map <- data.frame(
      rID = rownames(temp_geno),
      scanID = 1:length(rownames(temp_geno))
    )

    # Now reasign names
    rownames(temp_geno) <- temp_line_map$scanID
    temp_adj_means <- temp_line_map |>
      dplyr::left_join(temp_adj_means, by = "rID") |>
      tidyr::drop_na(scanID)

    # Subset marker_map
    temp_marker_map <- marker_map |>
      dplyr::filter(snp_id %in% colnames(temp_geno))

    # Make geontype data object for GENISIS
    geno_genisis <- GWASTools::MatrixGenotypeReader(
      genotype = t(temp_geno),
      snpID = temp_marker_map$variant.id,
      chromosome = temp_marker_map$chr,
      position = temp_marker_map$pos,
      scanID = temp_line_map$scanID
    )

    # Make a genotype data object
    geno_data_genisis <- GWASTools::GenotypeData(geno_genisis)

    # Calculate kinship matrix
    temp_kinship_file <- sommer::A.mat(temp_geno - 1)

    # Do PCA
    temp_pcs_geno <- prcomp(temp_geno)

    # Get list of PCs to include
    temp_pcs_keep <- summary(temp_pcs_geno)$importance |>
      t() |>
      as.data.frame() |>
      dplyr::filter(`Proportion of Variance` >= 0.03)
    temp_pcs_keep <- rownames(temp_pcs_keep)

    # Make a new covar file for the final scanannot file
    covar_file_scanAnnot <- temp_pcs_geno$x |>
      as.data.frame() |>
      tibble::rownames_to_column(var = "scanID")

    # Pull first x PCs
    covar_file_scanAnnot <- as.data.frame(covar_file_scanAnnot[, c("scanID", temp_pcs_keep)])
    covar_file_scanAnnot$scanID <- as.integer(covar_file_scanAnnot$scanID)

    # Make scanannot
    scanAnnot <- temp_adj_means |>
      dplyr::left_join(covar_file_scanAnnot, by = "scanID") |>
      dplyr::select(scanID, all_of(temp_pcs_keep), predicted.value) |>
      GWASTools::ScanAnnotationDataFrame()

    # Temp Null model
    temp_null <- try(
      GENESIS::fitNullModel(
        scanAnnot,
        outcome = "predicted.value",
        covars = temp_pcs_keep,
        cov.mat = temp_kinship_file
      )
    )

    # If try fail
    if (class(temp_null) == "try-error") {
      # Say something
      print(paste("GWAS of", j, "in", i, "failed. Skipping..."))

      # Remove
      remove(
        covar_file_scanAnnot,
        geno_data_genisis,
        scanAnnot,
        temp_adj_means,
        temp_geno,
        temp_kinship_file,
        temp_line_map,
        temp_marker_map,
        temp_null,
        temp_pcs_geno,
        temp_pheno_sub,
        temp_pcs_keep
      )

      # Skip
      next()
    }

    # Temp iterator
    temp_iterator <- GWASTools::GenotypeBlockIterator(geno_data_genisis, snpBlock = 1000)

    # Temp GWAS model
    temp_gwas <- suppressWarnings(GENESIS::assocTestSingle(temp_iterator, null.model = temp_null))

    # Temp output
    temp_out <- temp_marker_map |>
      dplyr::mutate(chr = as.character(chr)) |>
      dplyr::left_join(temp_gwas, by = c("chr", "variant.id", "pos")) |>
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
        break_bft = ifelse(
          neg_log_pval >= -log10(0.05 / nrow(temp_gwas)),
          1,
          0
        )
      )

    # temp plot
    temp_plot <- plot_gwas(
      temp_out,
      title = paste0(j, " in ", i, " (n = ", nrow(temp_pheno_sub), "; nPC = ", length(temp_pcs_keep), ")"),
      threshold = -log10(0.05 / nrow(temp_gwas))
    )

    # Save image
    ggplot2::ggsave(
      filename = file.path(
        temp_filepath,
        paste(j, "_in_", gsub(":", "_", i), ".jpg", sep = "")
      ),
      plot = temp_plot,
      width = 11,
      height = 8.5,
      units = "in",
      dpi = 180
    )

    # Save results
    write.csv(
      temp_out,
      file.path(
        temp_filepath,
        paste(j, "_in_", gsub(":", "_", i), ".csv", sep = "")
      )
    )

    # Make message for system time
    print_message <- paste("### System Time =", Sys.time(), "###")
    border <- paste(rep("#", nchar(print_message)), collapse = "")
    print(border)
    print(print_message)
    print(border)

    # Remove temporary stuff
    remove(
      covar_file_scanAnnot,
      geno_data_genisis,
      geno_genisis,
      scanAnnot,
      temp_adj_means,
      temp_geno,
      temp_gwas,
      temp_iterator,
      temp_kinship_file,
      temp_line_map,
      temp_null,
      temp_out,
      temp_pcs_geno,
      temp_pheno_sub,
      temp_pcs_keep,
      temp_plot,
      temp_marker_map
    )
  }

  # Clean up
  remove(
    temp_pheno,
    temp_filepath
  )
}
