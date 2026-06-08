# Read in RDS
PreDCS <- readRDS("PreDCS_GWAS/pca_results.RDS")
PostDCS <- readRDS("PostDCS_GWAS/pca_results.RDS")

# Format Base PCA Data
PreDCS_PCA <- PreDCS$pcs_geno |>
  tibble::rownames_to_column(var = "scanID") |>
  as.data.frame() |>
  dplyr::mutate(scanID = as.numeric(scanID)) |>
  dplyr::left_join(PreDCS$line_map, by = "scanID") |>
  dplyr::select(rID, dplyr::everything())

PostDCS_PCA <- PostDCS$pcs_geno |>
  tibble::rownames_to_column(var = "scanID") |>
  as.data.frame() |>
  dplyr::mutate(scanID = as.numeric(scanID)) |>
  dplyr::left_join(PostDCS$line_map, by = "scanID") |>
  dplyr::select(rID, dplyr::everything())

# Get keyfile
keyfile_path <- file.path(getwd(), "../Raw_Data/virus_mapping_2025.tsv")

# Define all possible classes
all_shapes <- c(
  "Standard",
  "CoAxium",
  "Clearfield",
  "Sawfly",
  "Doubled Haploid",
  "CoAxium, Clearfield",
  "CoAxium, Sawfly",
  "CoAxium, Doubled Haploid",
  "Clearfield, Sawfly",
  "Clearfield, Doubled Haploid",
  "Sawfly, Doubled Haploid",
  "CoAxium, Clearfield, Sawfly",
  "CoAxium, Clearfield, Doubled Haploid",
  "CoAxium, Sawfly, Doubled Haploid",
  "Clearfield, Sawfly, Doubled Haploid",
  "CoAxium, Clearfield, Sawfly, Doubled Haploid"
)

# Define all possible shapes
shape_palette <- c(
  "Standard" = 16,
  "CoAxium" = 15,
  "Clearfield" = 17,
  "Sawfly" = 4,
  "Doubled Haploid" = 3,
  "CoAxium, Clearfield" = 7,
  "CoAxium, Sawfly" = 8,
  "CoAxium, Doubled Haploid" = 9,
  "Clearfield, Sawfly" = 10,
  "Clearfield, Doubled Haploid" = 12,
  "Sawfly, Doubled Haploid" = 13,
  "CoAxium, Clearfield, Sawfly" = 14,
  "CoAxium, Clearfield, Doubled Haploid" = 18,
  "CoAxium, Sawfly, Doubled Haploid" = 5,
  "Clearfield, Sawfly, Doubled Haploid" = 6,
  "CoAxium, Clearfield, Sawfly, Doubled Haploid" = 11
)

# Define base classifications
base_class_levels <- c(
  "Red Winter",
  "White Winter",
  "Mixed Winter",
  "Unclassified / Single Trait"
)

# Define the colorblind-friendly color palette
color_palette <- c(
  "Red Winter" = "#D55E00", # Vermilion
  "White Winter" = "#56B4E9", # Sky Blue
  "Mixed Winter" = "#009E73", # Bluish Green
  "Unclassified / Single Trait" = "#000000" # Black
)

# Define base names
base_names <- c(
  "Red Winter",
  "White Winter",
  "Mixed Winter"
)

# Unclassified lines that have only trait designations
traits_only <- all_shapes[all_shapes != "Standard"]

# Define all possible classes
all_market_classes <- c(
  "Unclassified",
  base_names,
  traits_only,
  paste0(
    as.vector(
      outer(base_names, traits_only, paste, sep = " (")
    ),
    ")"
  )
)

# Create function to make classifications
format_wheat_classes <- function(pca_data) {
  pca_data |>
    # Mutate new classifications
    dplyr::mutate(
      # Find out if its a CO line
      `Experimental Line` = startsWith(FullSampleName, "CO"),
      # Pull all leters except CO
      Clean_Traits = ifelse(`Experimental Line` == TRUE, gsub("^CO|[0-9]+", "", FullSampleName), NA_character_),
      # Define the base class
      Base_Class = dplyr::case_when(
        `Experimental Line` == TRUE & grepl("R", Clean_Traits) ~ "Red Winter",
        `Experimental Line` == TRUE & grepl("W", Clean_Traits) ~ "White Winter",
        `Experimental Line` == TRUE & grepl("M", Clean_Traits) ~ "Mixed Winter",
        .default = ""
      ),
      # Now paste togeter subclasses
      Stacked_Traits = paste0(
        ifelse(`Experimental Line` == TRUE & grepl("A", Clean_Traits), "CoAxium, ", ""),
        ifelse(`Experimental Line` == TRUE & grepl("C", Clean_Traits), "Clearfield, ", ""),
        ifelse(`Experimental Line` == TRUE & grepl("SF", Clean_Traits), "Sawfly, ", ""),
        ifelse(`Experimental Line` == TRUE & grepl("D", Clean_Traits), "Doubled Haploid, ", "")
      ),
      # Replace the ,  if it is in the stacked trait name
      Stacked_Traits = sub(", $", "", Stacked_Traits),
      # Now add market class
      raw_Market_Class = dplyr::case_when(
        `Experimental Line` == FALSE ~ "Unclassified",
        Base_Class != "" & Stacked_Traits != "" ~ paste0(Base_Class, " (", Stacked_Traits, ")"),
        Base_Class != "" & Stacked_Traits == "" ~ Base_Class,
        Base_Class == "" & Stacked_Traits != "" ~ Stacked_Traits,
        .default = NA_character_
      ),
      # Add the group color
      raw_Color_Group = dplyr::case_when(
        `Experimental Line` == FALSE ~ "Unclassified / Single Trait",
        Base_Class == "Red Winter" ~ "Red Winter",
        Base_Class == "White Winter" ~ "White Winter",
        Base_Class == "Mixed Winter" ~ "Mixed Winter",
        `Experimental Line` == TRUE & Base_Class == "" ~ "Unclassified / Single Trait",
        .default = NA_character_
      ),
      # Add the shape associated with the class
      raw_Shape = dplyr::if_else(Stacked_Traits == "" | is.na(Stacked_Traits), "Standard", Stacked_Traits),
      # Make market class a factor
      `Market Class` = factor(raw_Market_Class, levels = all_market_classes),
      # Make plot color group
      Plot_Color_Group = factor(raw_Color_Group, levels = base_class_levels),
      # Make a shape
      Plot_Shape_Group = factor(raw_Shape, levels = all_shapes)
    ) |>
    # Select all those different columns we just made and put them to the front of the object
    dplyr::select(
      `Experimental Line`,
      `Market Class`,
      Plot_Color_Group,
      Plot_Shape_Group,
      dplyr::everything(),
      -Clean_Traits,
      -Base_Class,
      -Stacked_Traits,
      -raw_Market_Class,
      -raw_Color_Group,
      -raw_Shape
    )
}

# Now check for keyfile
if (file.exists(keyfile_path)) {
  # Read in keyfile
  keyfile <- read.table(keyfile_path, header = TRUE) |>
    dplyr::select(RealFullSampleName, FullSampleName) |>
    dplyr::rename(rID = FullSampleName, FullSampleName = RealFullSampleName) |>
    dplyr::distinct()

  # Format PreDCS PCA
  PreDCS_PCA <- PreDCS_PCA |>
    dplyr::left_join(keyfile, by = "rID") |>
    dplyr::select(FullSampleName, dplyr::everything()) |>
    format_wheat_classes()

  # Format PostDCS PCA
  PostDCS_PCA <- PostDCS_PCA |>
    dplyr::left_join(keyfile, by = "rID") |>
    dplyr::select(FullSampleName, dplyr::everything()) |>
    format_wheat_classes()

  # Library packages for visualization
  library(ggplot2)
  library(patchwork)

  # Calculate percent variance for plots
  pre_vars <- round(PreDCS$pve_geno[2, c("PC1", "PC2", "PC3")] * 100, 1)
  post_vars <- round(PostDCS$pve_geno[2, c("PC1", "PC2", "PC3")] * 100, 1)

  # Make plot helper function
  plot_pca <- function(
    df, # PCA dataframe with columns
    x_col, # What is the X
    y_col, # What is the Y
    x_var, # What is the x percent variance
    y_var, # what is the y percent variance
    title, # What is the title
    tag = NULL, # What is the tag level
    show_legend = TRUE, # Do you want to see the legend in the image
    color_pal = color_palette, # Allow custom palette overrides
    shape_pal = shape_palette # Allow custom shape overrides
  ) {
    # Render plot
    p <- ggplot(
      df,
      aes(
        x = .data[[x_col]],
        y = .data[[y_col]],
        color = Plot_Color_Group,
        shape = Plot_Shape_Group
      )
    ) +
      geom_point(size = 1, alpha = 0.7) +
      scale_color_manual(
        values = color_pal,
        name = "Market Class",
        guide = guide_legend(override.aes = list(shape = 16))
      ) +
      scale_shape_manual(
        values = shape_pal
      ) +
      labs(
        title = title,
        x = paste0(x_col, " (", x_var, "%)"),
        y = paste0(y_col, " (", y_var, "%)"),
        shape = "Subclassification",
        tag = tag
      ) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 11))

    # Suppress legend
    if (!show_legend) {
      p <- p + theme(legend.position = "none")
    }

    # Return plot
    return(p)
  }

  # Render PreDCS plots
  p_pre_12 <- plot_pca(PreDCS_PCA, "PC1", "PC2", pre_vars["PC1"], pre_vars["PC2"], "PC1 vs PC2", tag = "A.", show_legend = FALSE)
  p_pre_13 <- plot_pca(PreDCS_PCA, "PC1", "PC3", pre_vars["PC1"], pre_vars["PC3"], "PC1 vs PC3", show_legend = FALSE)
  p_pre_23 <- plot_pca(PreDCS_PCA, "PC2", "PC3", pre_vars["PC2"], pre_vars["PC3"], "PC2 vs PC3", show_legend = FALSE)

  # Render PostDCS plots
  p_post_12 <- plot_pca(PostDCS_PCA, "PC1", "PC2", post_vars["PC1"], post_vars["PC2"], "PC1 vs PC2", tag = "B.")
  p_post_13 <- plot_pca(PostDCS_PCA, "PC1", "PC3", post_vars["PC1"], post_vars["PC3"], "PC1 vs PC3")
  p_post_23 <- plot_pca(PostDCS_PCA, "PC2", "PC3", post_vars["PC2"], post_vars["PC3"], "PC2 vs PC3")

  # Combine plots
  final_plot <- (p_pre_12 | p_pre_13 | p_pre_23) /
    (p_post_12 | p_post_13 | p_post_23) +
    plot_layout(guides = "collect")

  # Save plot
  ggsave(
    filename = file.path(getwd(), "PCA_image/pre_and_post_DCS.jpg"),
    plot = final_plot,
    width = 10,
    height = 5,
    dpi = 320,
    units = "in",
    create.dir = TRUE
  )

  # Define colorblind-friendly simplified palettes
  simp_color_palette <- c(
    "Red" = "#D55E00", # Vermilion
    "White" = "#56B4E9", # Sky Blue
    "Mixed" = "#009E73", # Bluish Green
    "Unclassified" = "#000000" # Black
  )

  simp_shape_palette <- c(
    "Standard" = 16,
    "Doubled Haploid" = 3
  )

  # Create a function to override/simplify the data frames
  simplify_pca_classes <- function(df) {
    df |>
      dplyr::mutate(
        # Rename color groups
        Plot_Color_Group = dplyr::case_when(
          Plot_Color_Group == "Red Winter" ~ "Red",
          Plot_Color_Group == "White Winter" ~ "White",
          Plot_Color_Group == "Mixed Winter" ~ "Mixed",
          Plot_Color_Group == "Unclassified / Single Trait" ~ "Unclassified",
          .default = as.character(Plot_Color_Group)
        ),
        Plot_Color_Group = factor(Plot_Color_Group, levels = names(simp_color_palette)),

        # Collapse shapes into just DH and Standard
        Plot_Shape_Group = dplyr::if_else(
          grepl("Doubled Haploid", as.character(Plot_Shape_Group)),
          "Doubled Haploid",
          "Standard"
        ),
        Plot_Shape_Group = factor(Plot_Shape_Group, levels = names(simp_shape_palette))
      )
  }

  # Apply simplified classifications
  PreDCS_PCA_simp <- simplify_pca_classes(PreDCS_PCA)
  PostDCS_PCA_simp <- simplify_pca_classes(PostDCS_PCA)

  # Render Simplified PreDCS plots
  p_pre_12_simp <- plot_pca(PreDCS_PCA_simp, "PC1", "PC2", pre_vars["PC1"], pre_vars["PC2"], "PC1 vs PC2", tag = "A.", show_legend = FALSE, color_pal = simp_color_palette, shape_pal = simp_shape_palette)
  p_pre_13_simp <- plot_pca(PreDCS_PCA_simp, "PC1", "PC3", pre_vars["PC1"], pre_vars["PC3"], "PC1 vs PC3", show_legend = FALSE, color_pal = simp_color_palette, shape_pal = simp_shape_palette)
  p_pre_23_simp <- plot_pca(PreDCS_PCA_simp, "PC2", "PC3", pre_vars["PC2"], pre_vars["PC3"], "PC2 vs PC3", show_legend = FALSE, color_pal = simp_color_palette, shape_pal = simp_shape_palette)

  # Render Simplified PostDCS plots
  p_post_12_simp <- plot_pca(PostDCS_PCA_simp, "PC1", "PC2", post_vars["PC1"], post_vars["PC2"], "PC1 vs PC2", tag = "B.", color_pal = simp_color_palette, shape_pal = simp_shape_palette)
  p_post_13_simp <- plot_pca(PostDCS_PCA_simp, "PC1", "PC3", post_vars["PC1"], post_vars["PC3"], "PC1 vs PC3", color_pal = simp_color_palette, shape_pal = simp_shape_palette)
  p_post_23_simp <- plot_pca(PostDCS_PCA_simp, "PC2", "PC3", post_vars["PC2"], post_vars["PC3"], "PC2 vs PC3", color_pal = simp_color_palette, shape_pal = simp_shape_palette)

  # Combine simplified plots
  final_plot_simp <- (p_pre_12_simp | p_pre_13_simp | p_pre_23_simp) /
    (p_post_12_simp | p_post_13_simp | p_post_23_simp) +
    plot_layout(guides = "collect")

  # Save simplified plot
  ggsave(
    filename = file.path(getwd(), "PCA_image/pre_and_post_DCS_simplified.jpg"),
    plot = final_plot_simp,
    width = 10,
    height = 5,
    dpi = 320,
    units = "in",
    create.dir = TRUE
  )
}
