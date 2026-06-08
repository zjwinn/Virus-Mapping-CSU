# Read in PreDCS, PostDCS, and PostDCS Binomeal GWAS
predcs_gwas <- read.csv("PreDCS_GWAS/gwas_results_continuous_pre_yd.csv")
postdcs_gwas_c <- read.csv("PostDCS_GWAS/gwas_results_continuous_post_yd.csv")
postdcs_gwas_b <- read.csv("PostDCS_GWAS_binomial/gwas_results_binomial.csv")

# Add new classifires to new data
gwas <- predcs_gwas |>
  dplyr::mutate(TYPE = "PreDCS") |>
  dplyr::bind_rows(postdcs_gwas_c |> dplyr::mutate(TYPE = "PostDCS_C")) |>
  dplyr::bind_rows(postdcs_gwas_b |> dplyr::mutate(TYPE = "PostDCS_B"))

# Make gwas function
plot_gwas <- function(df,
                      chr = "chr_wheat",
                      pos = "pos",
                      y = "neg_log_pval",
                      threshold = NULL,
                      title = NULL,
                      lim = c(0, ceiling(max(df$neg_log_pval)))) {
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
    ggplot2::scale_y_continuous(limits = lim) +
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

# Make a function for QQ-plots and lamda
gwas_qq <- function(assoc, pval_col = "Score.pval") {
  # Get pval
  pval <- assoc[[pval_col]]

  # Keep only pvals that are complete
  pval <- pval[!is.na(pval) & pval > 0]

  # calculate lambda from chi-square quantiles
  chisq <- qchisq(pval, df = 1, lower.tail = FALSE)
  lambda <- median(chisq) / qchisq(0.5, df = 1)

  # Make dataframe
  n <- length(pval)
  x <- 1:n
  dat <- data.frame(
    obs   = sort(pval),
    exp   = x / (n + 1),
    upper = qbeta(0.975, x, n - x + 1),
    lower = qbeta(0.025, x, n - x + 1)
  )

  # Plot QQ plot
  qq <- ggplot2::ggplot(dat, ggplot2::aes(-log10(exp), -log10(obs))) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = -log10(upper), ymax = -log10(lower)),
      fill = "grey80"
    ) +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red") +
    ggplot2::annotate(
      "text",
      x = 0.5, y = max(-log10(dat$obs)), hjust = 0,
      label = sprintf("lambda[GC] == %.3f", lambda), parse = TRUE
    ) +
    ggplot2::xlab(expression(paste(-log[10], "(expected P)"))) +
    ggplot2::ylab(expression(paste(-log[10], "(observed P)"))) +
    ggplot2::theme_bw()

  # Return List
  list(lambda = lambda, qqplot = qq)
}

# Get lim
lim <- c(
  0,
  gwas |>
    dplyr::filter(TYPE == "PreDCS") |>
    dplyr::pull(neg_log_pval) |>
    max(na.rm = TRUE) |>
    ceiling()
)

# Plot first one
p1 <- plot_gwas(
  gwas |>
    dplyr::filter(trait == "GY" & TYPE == "PreDCS"),
  threshold = -log10(
    0.05 / nrow(
      gwas |>
        dplyr::filter(trait == "GY" & TYPE == "PreDCS")
    )
  ),
  title = "Grain Yield",
  lim = lim
) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank()
  )

# Plot second one
p2 <- plot_gwas(
  gwas |>
    dplyr::filter(trait == "TW" & TYPE == "PreDCS"),
  threshold = -log10(
    0.05 / nrow(
      gwas |>
        dplyr::filter(trait == "TW" & TYPE == "PreDCS")
    )
  ),
  title = "Test Weight",
  lim = lim
) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank()
  )

# Plot third one
p3 <- plot_gwas(
  gwas |>
    dplyr::filter(trait == "WSMV" & TYPE == "PreDCS"),
  threshold = -log10(
    0.05 / nrow(
      gwas |>
        dplyr::filter(trait == "WSMV" & TYPE == "PreDCS")
    )
  ),
  title = "WSMV Ratings (Transformed)",
  lim = lim
)

# Library
library(patchwork)

# Visualize
p4 <- p1 / p2 / p3

# Save
ggplot2::ggsave(
  filename = file.path(
    getwd(),
    "GWAS_visualizations",
    "PreDCS_GWAS_all_traits.jpg"
  ),
  plot = p4,
  width = 8,
  height = 8,
  units = "in",
  dpi = 320,
  create.dir = TRUE
)

# Now do qqplot and lambda
gwas_qq(gwas |> dplyr::filter(trait == "GY" & TYPE == "PreDCS"))
gwas_qq(gwas |> dplyr::filter(trait == "TW" & TYPE == "PreDCS"))
gwas_qq(gwas |> dplyr::filter(trait == "WSMV" & TYPE == "PreDCS"))

# Define SNPs to keep per trait
keep_snps <- predcs_gwas |>
  dplyr::filter(break_bft == 1) |>
  dplyr::group_by(trait, chr_wheat) |>
  dplyr::summarise(
    Proximal = snp_id[which.min(pos)],
    Distal = snp_id[which.max(pos)],
    Peak = snp_id[which.max(neg_log_pval)],
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(
    cols = c(Proximal, Peak, Distal),
    names_to = "type",
    values_to = "snp_id"
  ) |>
  dplyr::select(trait, snp_id, type)

# Filter predcs_gwas using trait + snp_id
summary_table <- predcs_gwas |>
  dplyr::filter(break_bft == 1) |>
  dplyr::inner_join(
    keep_snps,
    by = c("trait", "snp_id")
  ) |>
  dplyr::select(
    trait,
    chr_wheat,
    type,
    snp_id,
    pos,
    freq,
    neg_log_pval,
    Est,
    Est.SE,
    PVE
  ) |>
  dplyr::mutate(
    pos = (pos / 1000000),
    Est = ifelse(freq > 0.5, -Est, Est),
    freq = ifelse(freq > 0.5, 1 - freq, freq)
  ) |>
  dplyr::rename(
    Trait = trait,
    Chromosome = chr_wheat,
    `Relative Position` = type,
    SNP = snp_id,
    `Position (Mbp)` = pos,
    Frequency = freq,
    Effect = Est,
    SE = Est.SE,
    `Percent Variance` = PVE
  ) |>
  dplyr::mutate(
    Trait = dplyr::case_when(
      Trait == "GY" ~ "Grain Yield (Kg/Ha)",
      Trait == "TW" ~ "Test Weight (Kg/Hl)",
      Trait == "WSMV" ~ "WSMV Ratings (Transformed)"
    )
  )

# Write csv
write.csv(
  summary_table,
  file.path(getwd(), "GWAS_visualizations/PreDCS_gwas_table.csv"),
  row.names = FALSE
)

# Get lim
lim <- c(
  0,
  gwas |>
    dplyr::filter(TYPE == "PostDCS_C") |>
    dplyr::pull(neg_log_pval) |>
    max(na.rm = TRUE) |>
    ceiling()
)

# Plot first one
p1 <- plot_gwas(
  gwas |>
    dplyr::filter(trait == "GY" & TYPE == "PostDCS_C"),
  threshold = -log10(
    0.05 / nrow(
      gwas |>
        dplyr::filter(trait == "GY" & TYPE == "PostDCS_C")
    )
  ),
  title = "Grain Yield",
  lim = lim
) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank()
  )

# Plot second one
p2 <- plot_gwas(
  gwas |>
    dplyr::filter(trait == "TW" & TYPE == "PostDCS_C"),
  threshold = -log10(
    0.05 / nrow(
      gwas |>
        dplyr::filter(trait == "TW" & TYPE == "PostDCS_C")
    )
  ),
  title = "Test Weight",
  lim = lim
) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank()
  )

# Plot third one
p3 <- plot_gwas(
  gwas |>
    dplyr::filter(trait == "Virus" & TYPE == "PostDCS_C"),
  threshold = -log10(
    0.05 / nrow(
      gwas |>
        dplyr::filter(trait == "Virus" & TYPE == "PostDCS_C")
    )
  ),
  title = "Virus Ratings (Transformed)",
  lim = lim
)


# Visualize
p4 <- p1 / p2 / p3

# Save
ggplot2::ggsave(
  filename = file.path(
    getwd(),
    "GWAS_visualizations",
    "PostDCS_C_GWAS_all_traits.jpg"
  ),
  plot = p4,
  width = 8,
  height = 8,
  units = "in",
  dpi = 320,
  create.dir = TRUE
)

# Plot QQ
gwas_qq(gwas |> dplyr::filter(trait == "GY" & TYPE == "PostDCS_C"))
gwas_qq(gwas |> dplyr::filter(trait == "TW" & TYPE == "PostDCS_C"))
gwas_qq(gwas |> dplyr::filter(trait == "Virus" & TYPE == "PostDCS_C"))

# Define SNPs to keep per trait
keep_snps <- postdcs_gwas_c |>
  dplyr::filter(break_bft == 1) |>
  dplyr::group_by(trait, chr_wheat) |>
  dplyr::summarise(
    Proximal = snp_id[which.min(pos)],
    Distal = snp_id[which.max(pos)],
    Peak = snp_id[which.max(neg_log_pval)],
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(
    cols = c(Proximal, Peak, Distal),
    names_to = "type",
    values_to = "snp_id"
  ) |>
  dplyr::select(trait, snp_id, type)

# Filter predcs_gwas using trait + snp_id
summary_table <- postdcs_gwas_c |>
  dplyr::filter(break_bft == 1) |>
  dplyr::inner_join(
    keep_snps,
    by = c("trait", "snp_id")
  ) |>
  dplyr::select(
    trait,
    chr_wheat,
    type,
    snp_id,
    pos,
    freq,
    neg_log_pval,
    Est,
    Est.SE,
    PVE
  ) |>
  dplyr::mutate(
    pos = (pos / 1000000),
    Est = ifelse(freq > 0.5, -Est, Est),
    freq = ifelse(freq > 0.5, 1 - freq, freq)
  ) |>
  dplyr::rename(
    Trait = trait,
    Chromosome = chr_wheat,
    `Relative Position` = type,
    SNP = snp_id,
    `Position (Mbp)` = pos,
    Frequency = freq,
    Effect = Est,
    SE = Est.SE,
    `Percent Variance` = PVE
  ) |>
  dplyr::mutate(
    Trait = dplyr::case_when(
      Trait == "GY" ~ "Grain Yield (Kg/Ha)",
      Trait == "TW" ~ "Test Weight (Kg/Hl)",
      Trait == "Virus" ~ "Virus Ratings (Transformed)"
    )
  )

# Curate table
summary_table <- summary_table[
  c(2, 5, 8, 10:12, 14, 17),
]

# Write csv
write.csv(
  summary_table,
  file.path(getwd(), "GWAS_visualizations/PostDCS_C_gwas_table.csv"),
  row.names = FALSE
)

# Get lim
lim <- c(
  0,
  gwas |>
    dplyr::filter(TYPE == "PostDCS_B") |>
    dplyr::pull(neg_log_pval) |>
    max(na.rm = TRUE) |>
    ceiling()
)

# Plot first one
p1 <- plot_gwas(
  gwas |>
    dplyr::filter(trait == "Virus" & TYPE == "PostDCS_B"),
  threshold = -log10(
    0.05 / nrow(
      gwas |>
        dplyr::filter(trait == "Virus" & TYPE == "PostDCS_B")
    )
  ),
  title = "Virus Ratings (Binomial)",
  lim = lim
)

# Save image
ggplot2::ggsave(
  filename = file.path(
    getwd(),
    "GWAS_visualizations",
    "PostDCS_B_GWAS_all_traits.jpg"
  ),
  plot = p1,
  width = 8,
  height = 2.67,
  units = "in",
  dpi = 320,
  create.dir = TRUE
)

# see QQ plot
gwas_qq(gwas |> dplyr::filter(trait == "Virus" & TYPE == "PostDCS_B"))

# Define SNPs to keep per trait
keep_snps <- postdcs_gwas_b |>
  dplyr::filter(break_bft == 1) |>
  dplyr::group_by(trait, chr_wheat) |>
  dplyr::summarise(
    Proximal = snp_id[which.min(pos)],
    Distal = snp_id[which.max(pos)],
    Peak = snp_id[which.max(neg_log_pval)],
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(
    cols = c(Proximal, Peak, Distal),
    names_to = "type",
    values_to = "snp_id"
  ) |>
  dplyr::select(trait, snp_id, type)

# Filter predcs_gwas using trait + snp_id
summary_table <- postdcs_gwas_b |>
  dplyr::filter(break_bft == 1) |>
  dplyr::inner_join(
    keep_snps,
    by = c("trait", "snp_id")
  ) |>
  dplyr::select(
    chr_wheat,
    type,
    snp_id,
    pos,
    freq,
    neg_log_pval,
    Est,
    Est.SE,
    PVE
  ) |>
  dplyr::mutate(
    pos = (pos / 1000000),
    Est = ifelse(freq > 0.5, -Est, Est),
    freq = ifelse(freq > 0.5, 1 - freq, freq)
  ) |>
  dplyr::rename(
    Chromosome = chr_wheat,
    `Relative Position` = type,
    SNP = snp_id,
    `Position (Mbp)` = pos,
    Frequency = freq,
    Effect = Est,
    SE = Est.SE,
    `Percent Variance` = PVE
  )

# Curate table
summary_table <- summary_table[
  c(2, 4:6, 8, 11),
]

# Write csv
write.csv(
  summary_table,
  file.path(getwd(), "GWAS_visualizations/PostDCS_B_gwas_table.csv"),
  row.names = FALSE
)
