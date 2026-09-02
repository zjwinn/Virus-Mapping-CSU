# Initialize required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(lubridate)
library(scales)
library(sf)
library(tigris)
library(patchwork)
library(svglite)

# Enable tigris geometry caching
options(tigris_use_cache = TRUE)

# Save original working directory and ensure safe restoration
orig_dir <- getwd()
on.exit(setwd(orig_dir))

# Set working directory
setwd(file.path(orig_dir, "virus_frequency_analysis/"))

# Ingest and format diagnostic data
raw_data <- read_excel(
  "Diagnostic Results_Mason.xlsx",
  sheet = "Table_all diagnostics",
  skip = 3
) |>
  select(`Date Received`, Location, WSMV, TriMV) |>
  mutate(
    Year = year(`Date Received`),
    State = case_when(
      WSMV == "-" & TriMV == "-" ~ "Uninfected",
      WSMV == "+" & TriMV == "-" ~ "WSMV Only",
      WSMV == "-" & TriMV == "+" ~ "TriMV Only",
      WSMV == "+" & TriMV == "+" ~ "Co-infected",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(State) & !is.na(Year)) |>
  mutate(
    State = factor(
      State,
      levels = c("Uninfected", "WSMV Only", "TriMV Only", "Co-infected")
    )
  ) |>
  as.data.frame()

# Calculate empirical frequencies and zero-fill missing factors
empirical_df <- raw_data |>
  group_by(Year, State) |>
  summarise(Count = n(), .groups = "drop") |>
  complete(Year, State, fill = list(Count = 0)) |>
  group_by(Year) |>
  mutate(
    Total_Samples = sum(Count),
    Proportion = Count / Total_Samples,
    Label = ifelse(Proportion > 0, percent(Proportion, accuracy = 1), "")
  ) |>
  ungroup()

# Isolate sample sizes for x-axis labeling
axis_labels <- empirical_df |>
  select(Year, Total_Samples) |>
  distinct() |>
  arrange(Year)

# Parse and aggregate spatial location data
spatial_df <- raw_data |>
  filter(!is.na(Location) & tolower(Location) != "unknown") |>
  mutate(
    Clean_Location = str_remove(Location, " County"),
    County_Name = trimws(str_split_i(Clean_Location, ",", 1)),
    State_Abbr = trimws(str_split_i(Clean_Location, ",", 2))
  ) |>
  filter(!is.na(State_Abbr)) |>
  group_by(Year, State_Abbr, County_Name) |>
  summarise(Total_Assessments = n(), .groups = "drop")

# Retrieve county and state geometries
target_states <- unique(spatial_df$State_Abbr)
county_shapes <- counties(state = target_states, cb = TRUE, class = "sf")
state_shapes <- states(cb = TRUE, class = "sf") |> filter(STUSPS %in% target_states)

# Join FIPS codes to spatial data
fips_dictionary <- fips_codes |>
  select(state, state_code) |>
  distinct()
spatial_df <- spatial_df |> left_join(fips_dictionary, by = c("State_Abbr" = "state"))

# Generate complete spatiotemporal grid and merge count data
temporal_spatial_grid <- crossing(
  Year = unique(spatial_df$Year),
  state_code = county_shapes$STATEFP,
  County_Name = county_shapes$NAME
)

completed_spatial_df <- temporal_spatial_grid |>
  left_join(spatial_df, by = c("Year", "state_code", "County_Name")) |>
  mutate(Total_Assessments = replace_na(Total_Assessments, 0))

map_data <- county_shapes |>
  left_join(completed_spatial_df, by = c("STATEFP" = "state_code", "NAME" = "County_Name"))

# Generate temporal visualization
pub_plot_1 <- ggplot(empirical_df, aes(x = Year, y = Proportion, fill = State)) +
  geom_col(position = "stack", width = 0.7, color = "black") +
  geom_text(
    aes(label = Label),
    position = position_stack(vjust = 0.5),
    size = 4,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("gray90", "#E69F00", "#56B4E9", "#D55E00")) +
  scale_x_continuous(
    breaks = axis_labels$Year,
    labels = paste0(axis_labels$Year, "\n(", axis_labels$Total_Samples, ")")
  ) +
  scale_y_continuous(
    labels = percent_format(),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Assessment Year",
    y = "Proportion of Total Samples",
    fill = "Infection State"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

# Generate faceted spatial visualization
pub_plot_2 <- ggplot() +
  geom_sf(data = map_data, aes(fill = Total_Assessments), color = "black", linewidth = 0.2) +
  geom_sf(data = state_shapes, fill = NA, color = "black", linewidth = 1.2) +
  facet_wrap(~Year, ncol = 2) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    name = "Assessments",
    na.value = "white"
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(margin = margin(b = 5))
  )

# Compile and export manuscript composite
pub_combined <- (pub_plot_1 + pub_plot_2) +
  plot_layout(widths = c(1, 1.5)) +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = ".",
    theme = theme(plot.tag = element_text(size = 18, face = "bold"))
  )

ggsave(
  "Manuscript_Wheat_Virus_Color.jpeg",
  plot = pub_combined,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300
)

# Generate temporal visualization
poster_plot_1 <- ggplot(empirical_df, aes(x = Year, y = Proportion, fill = State)) +
  geom_col(position = "stack", width = 0.7, color = "black") +
  geom_text(
    aes(label = Label, color = State),
    position = position_stack(vjust = 0.5),
    size = 5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("#E6E6E6", "#A6A6A6", "#595959", "#1A1A1A")) +
  scale_color_manual(values = c("black", "black", "white", "white")) +
  scale_x_continuous(
    breaks = axis_labels$Year,
    labels = paste0(axis_labels$Year, "\n(", axis_labels$Total_Samples, ")")
  ) +
  scale_y_continuous(
    labels = percent_format(),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Assessment Year",
    y = "Proportion of Total Samples",
    fill = "Infection State"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

# Generate faceted spatial visualization
poster_plot_2 <- ggplot() +
  geom_sf(data = map_data, aes(fill = Total_Assessments), color = "black", linewidth = 0.2) +
  geom_sf(data = state_shapes, fill = NA, color = "black", linewidth = 1.2) +
  facet_wrap(~Year, ncol = 2) +
  scale_fill_gradient(
    low = "gray95",
    high = "black",
    name = "Assessments",
    na.value = "white"
  ) +
  theme_void(base_size = 16) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(margin = margin(b = 5))
  )

# Compile and export poster composite
poster_combined <- (poster_plot_1 + poster_plot_2) +
  plot_layout(widths = c(1, 1.5)) +
  plot_annotation(
    title = "Wheat Curlmite Vectored Virus Frequencies are Shifting in Colorado",
    subtitle = "Data from Co-PI Roberts et al",
    tag_levels = "A",
    tag_suffix = ".",
    theme = theme(
      plot.title = element_text(size = 24),
      plot.subtitle = element_text(size = 18),
      plot.tag = element_text(size = 20)
    )
  )

ggsave(
  "Poster_Wheat_Virus_BW.svg",
  plot = poster_combined,
  width = 14,
  height = 8,
  units = "in",
  device = "svg"
)

ggsave(
  "Poster_Wheat_Virus_BW.jpeg",
  plot = poster_combined,
  width = 14,
  height = 8,
  units = "in",
  dpi = 300
)

# Set working directory
setwd(file.path(getwd(), ".."))
