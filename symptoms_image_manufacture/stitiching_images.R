library(magick)
library(cowplot)
library(ggplot2)

# Load and normalize raw images
img_A <- image_normalize(image_read("IMG_3737.tiff"))
img_B <- image_normalize(image_read("IMG_8246.tiff"))
img_C <- image_normalize(image_read("IMG_7658.tiff"))

# Center-crop an image to a specific aspect ratio
crop_to_ar <- function(img, target_ar) {
  info <- image_info(img)
  current_ar <- info$width / info$height
  if (current_ar > target_ar) {
    new_w <- round(info$height * target_ar)
    new_h <- info$height
  } else {
    new_w <- info$width
    new_h <- round(info$width / target_ar)
  }
  image_crop(img, paste0(new_w, "x", new_h), gravity = "center")
}

# Crop images to exactly fit an 11x8.5 grid layout (1:1.3 relative height)
img_A <- crop_to_ar(img_A, 5.5 / 3.7)
img_B <- crop_to_ar(img_B, 5.5 / 3.7)
img_C <- crop_to_ar(img_C, 11 / 4.8)

# Function to add a white label with a bold black drop-shadow for max legibility
add_better_label <- function(plot, lab) {
  plot + 
    draw_plot_label(lab, x = 0.023, y = 0.977, size = 24, colour = "black") +
    draw_plot_label(lab, x = 0.020, y = 0.980, size = 24, colour = "white")
}

# Convert magick images to ggplot objects and add the shadow labels
p_A <- add_better_label(ggdraw() + draw_image(img_A), "A.")
p_B <- add_better_label(ggdraw() + draw_image(img_B), "B.")
p_C <- add_better_label(ggdraw() + draw_image(img_C), "C.")

# Assemble the grid (labels are already baked into the individual plots now!)
top_row <- plot_grid(p_A, p_B, ncol = 2)
final_plot <- plot_grid(top_row, p_C, ncol = 1, rel_heights = c(1, 1.3))

# Save the plot at 300 dpi
ggsave("combined_stitched_image.jpg", final_plot, width = 11, height = 8.5, units = "in", dpi = 300)

# Compress the resulting 300 dpi jpg to reduce file size
image_write(image_read("combined_stitched_image.jpg"), "combined_stitched_image.jpg", quality = 75)
