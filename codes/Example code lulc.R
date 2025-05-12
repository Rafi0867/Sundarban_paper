#  ==========================    PREAAMBLE    ==========================
if(!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,
  ggplot2,
  terra,
  raster,
  mapview,
  dplyr,
  sf,
  lubridate,
  downloader,
  ncdf4
)



data <- rast("raster_data_location")
plot(data)

shapes <- vect("shape_files_location")
plot(shapes)



# Reproject shapefile if CRS doesn't match
if (!crs(data) == crs(shapes)) {
  shapes <- project(shapes, crs(data))
}

# Crop and mask
r_cropped <- crop(data, shapes)
r_masked <- mask(r_cropped, shapes)

plot(r_masked)


# Convert to data frame for ggplot
r_df <- as.data.frame(r_masked, xy = TRUE)
colnames(r_df)[3] <- "value"  # rename raster column



# Plot
ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = factor(value))) +
  geom_sf(data = st_as_sf(shapes), fill = NA, color = NA, size = 0.3) +
  scale_fill_manual(
    name = "Class",
    values = c("0" = "blue", "1" = "#008000", "2" = "orange", "3" = "yellow", "4" = "red"),
    labels = c("Water", "Mangrove", "Mixed", "Shrub", "Built-up")  # Edit if needed
  ) +
  coord_sf(expand = FALSE) +
  labs(
    title = "Land Cover Classification in Sundarbans",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )
