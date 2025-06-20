# install required packages----
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  terra, # handle raster data
  tidyterra, # handle and visualize raster data
  raster, # handle raster data
  exactextractr, # fast extractions
  sf, # vector data operations
  dplyr, # data wrangling
  tidyr, # data wrangling
  data.table, # data wrangling
  prism, # download PRISM data
  tictoc, # timing codes
  tigris, # to get county sf
  ggplot2,
  exactextractr,
  spatialEco,
  lulcc
)


# Load packages
library(terra)
library(readxl)
library(tidyverse)
library(lubridate)
library(patchwork)
library(viridis)
library(lulcc)
library(spatialEco)


# Define file paths for your LULC TIFFs
lulc_files <- list(
  "2000" = "Data/sundarban tifs/Sundarban_rf_00.tif",
  "2005" = "Data/sundarban tifs/Sundarban_rf_05.tif",
  "2010" = "Data/sundarban tifs/Sundarban_rf_10.tif",
  "2015" = "Data/sundarban tifs/Sundarban_rf_15.tif",
  "2020" = "Data/sundarban tifs/Sundarban_rf_21.tif"
)

# Define your LULC classes (adjust if your labels are different or non-sequential)
# 0=Water, 1=Forest, 2=Built-up, 3=Aquaculture, 4=Others
lulc_class_names <- c("Water", "Forest", "Built-up", "Aquaculture", "Others")
num_lulc_classes <- length(lulc_class_names) # Should be 5

# Output file prefixes for predictions
output_prefix <- "predicted_lulc_"

# --- Functions ---

#' Calculate Transition Matrix
#'
#' Calculates the transition matrix between two LULC raster layers.
#'
#' @param lulc_prev A raster object (SpatRaster or RasterLayer) of the previous LULC.
#' @param lulc_curr A raster object (SpatRaster or RasterLayer) of the current LULC.
#' @param num_classes The total number of LULC classes.
#' @return A matrix representing the transition probabilities.
#'         Rows are 'from' classes, columns are 'to' classes.
calculate_transition_matrix <- function(lulc_prev, lulc_curr, num_classes) {
  # Convert rasters to vectors, removing NA values if any
  vals_prev <- as.vector(lulc_prev)
  vals_curr <- as.vector(lulc_curr)
  
  # Ensure both vectors are of the same length and filter out NAs
  # Filter out NA values from both vectors simultaneously
  valid_indices <- which(!is.na(vals_prev) & !is.na(vals_curr))
  vals_prev_filtered <- vals_prev[valid_indices]
  vals_curr_filtered <- vals_curr[valid_indices]
  
  # Initialize transition count matrix
  transition_counts <- matrix(0, nrow = num_classes, ncol = num_classes)
  
  # Populate transition counts
  for (i in 1:length(vals_prev_filtered)) {
    from_class <- vals_prev_filtered[i]
    to_class <- vals_curr_filtered[i]
    # Ensure classes are within valid range (0 to num_classes-1)
    if (from_class >= 0 && from_class < num_classes && to_class >= 0 && to_class < num_classes) {
      transition_counts[from_class + 1, to_class + 1] <- transition_counts[from_class + 1, to_class + 1] + 1
    }
  }
  
  # Convert counts to probabilities
  transition_matrix <- matrix(0, nrow = num_classes, ncol = num_classes)
  for (i in 1:num_classes) {
    row_sum <- sum(transition_counts[i, ])
    if (row_sum > 0) {
      transition_matrix[i, ] <- transition_counts[i, ] / row_sum
    } else {
      # If a class never appeared, assume it stays the same (identity)
      transition_matrix[i, i] <- 1.0
    }
  }
  
  # Add row and column names for better readability
  rownames(transition_matrix) <- paste0("From_", lulc_class_names)
  colnames(transition_matrix) <- paste0("To_", lulc_class_names)
  
  return(transition_matrix)
}


#' Predict LULC using Markov Chain
#'
#' Predicts the next LULC map using a Markov chain transition matrix.
#'
#' @param current_lulc_raster A raster object (SpatRaster or RasterLayer) of the current LULC.
#' @param transition_matrix The 2D transition probability matrix (from_class, to_class).
#' @param num_classes The total number of LULC classes.
#' @return A new raster object with the predicted LULC classes.
predict_lulc_markov <- function(current_lulc_raster, transition_matrix, num_classes) {
  # Create an empty raster for the predicted output, with the same properties
  # as the input raster. Using terra's `rast` for efficiency.
  predicted_lulc_raster <- terra::rast(current_lulc_raster)
  
  # Get the values of the current LULC raster
  current_lulc_values <- terra::values(current_lulc_raster)
  
  # Initialize predicted values vector
  predicted_lulc_values <- numeric(length(current_lulc_values))
  
  # Iterate through each pixel
  for (i in seq_along(current_lulc_values)) {
    current_class <- current_lulc_values[i]
    
    if (!is.na(current_class) && current_class >= 0 && current_class < num_classes) {
      # Get probabilities for the current class (add 1 because R is 1-indexed)
      probabilities <- transition_matrix[current_class + 1, ]
      
      # Predict based on highest probability (deterministic)
      # For stochastic prediction, you would use sample() function here
      predicted_class <- which.max(probabilities) - 1 # Subtract 1 to get 0-indexed class
      
      predicted_lulc_values[i] <- predicted_class
    } else {
      # Keep NA values as NA, or assign a 'no data' value
      predicted_lulc_values[i] <- NA
    }
  }
  
  # Assign the predicted values back to the raster
  terra::values(predicted_lulc_raster) <- predicted_lulc_values
  
  return(predicted_lulc_raster)
}

# --- Main Script ---

# 1. Load LULC data
lulc_rasters <- list()
for (year_str in names(lulc_files)) {
  filepath <- lulc_files[[year_str]]
  cat(paste0("Loading LULC for ", year_str, " from ", filepath, "...\n"))
  # Use terra::rast for loading - generally faster and more modern
  lulc_rasters[[year_str]] <- terra::rast(filepath)
  # Ensure integer values if they are not already
  lulc_rasters[[year_str]] <- round(lulc_rasters[[year_str]])
}

# 2. Calculate average transition matrix from historical data
all_transition_matrices <- list()
years_sorted <- sort(as.numeric(names(lulc_rasters)))

for (i in 1:(length(years_sorted) - 1)) {
  prev_year <- as.character(years_sorted[i])
  curr_year <- as.character(years_sorted[i+1])
  
  cat(paste0("Calculating transition for ", prev_year, " to ", curr_year, "...\n"))
  
  # Ensure rasters are aligned and have the same extent/resolution
  # This is crucial. If they are not identical, resample the current to the previous.
  if (crs(lulc_rasters[[prev_year]]) != crs(lulc_rasters[[curr_year]])) {
    warning(paste0("CRSs mismatch for ", prev_year, " and ", curr_year, ". Re-projecting ", curr_year, " to match."))
    lulc_rasters[[curr_year]] <- terra::project(lulc_rasters[[curr_year]], lulc_rasters[[prev_year]])
  }
  if (!all(ext(lulc_rasters[[prev_year]]) == ext(lulc_rasters[[curr_year]]))) {
    warning(paste0("Extents mismatch for ", prev_year, " and ", curr_year, ". Cropping/extending ", curr_year, " to match."))
    lulc_rasters[[curr_year]] <- terra::extend(terra::crop(lulc_rasters[[curr_year]], lulc_rasters[[prev_year]]), lulc_rasters[[prev_year]])
  }
  if (!all(res(lulc_rasters[[prev_year]]) == res(lulc_rasters[[curr_year]]))) {
    warning(paste0("Resolutions mismatch for ", prev_year, " and ", curr_year, ". Resampling ", curr_year, " to match."))
    lulc_rasters[[curr_year]] <- terra::resample(lulc_rasters[[curr_year]], lulc_rasters[[prev_year]], method = "near")
  }
  
  matrix_i <- calculate_transition_matrix(
    lulc_prev = lulc_rasters[[prev_year]],
    lulc_curr = lulc_rasters[[curr_year]],
    num_lulc_classes
  )
  all_transition_matrices[[paste0(prev_year, "_", curr_year)]] <- matrix_i
}

# Average the transition matrices
# Convert list of matrices to 3D array for easier averaging
if (length(all_transition_matrices) > 0) {
  avg_transition_matrix <- Reduce("+", all_transition_matrices) / length(all_transition_matrices)
} else {
  stop("No historical transitions to average. Check your LULC files and years.")
}

cat("\nAverage Transition Matrix:\n")
print(avg_transition_matrix)

# 3. Predict LULC for 2025, 2030, and 2035

current_lulc_raster <- lulc_rasters[["2020"]]

# Predict 2025
cat("\nPredicting LULC for 2025...\n")
predicted_lulc_2025 <- predict_lulc_markov(current_lulc_raster, avg_transition_matrix, num_lulc_classes)
output_filename_2025 <- paste0(output_prefix, "2025.tif")
terra::writeRaster(predicted_lulc_2025, filename = output_filename_2025, overwrite = TRUE, datatype = "INT1U") # INT1U for 0-4 classes
cat(paste0("2025 LULC predicted and saved to ", output_filename_2025, "\n"))

# Predict 2030 (based on 2025 prediction)
cat("Predicting LULC for 2030...\n")
predicted_lulc_2030 <- predict_lulc_markov(predicted_lulc_2025, avg_transition_matrix, num_lulc_classes)
output_filename_2030 <- paste0(output_prefix, "2030.tif")
terra::writeRaster(predicted_lulc_2030, filename = output_filename_2030, overwrite = TRUE, datatype = "INT1U")
cat(paste0("2030 LULC predicted and saved to ", output_filename_2030, "\n"))

# Predict 2035 (based on 2030 prediction)
cat("Predicting LULC for 2035...\n")
predicted_lulc_2035 <- predict_lulc_markov(predicted_lulc_2030, avg_transition_matrix, num_lulc_classes)
output_filename_2035 <- paste0(output_prefix, "2035.tif")
terra::writeRaster(predicted_lulc_2035, filename = output_filename_2035, overwrite = TRUE, datatype = "INT1U")
cat(paste0("2035 LULC predicted and saved to ", output_filename_2035, "\n"))

cat("\nPrediction complete.\n")

# # --- Optional: Visualize a predicted map (requires a display) ---
plot(predicted_lulc_2035, main = "Predicted LULC 2035")
#You might want to assign colors for better visualization
colors <- c("skyblue", "darkgreen", "red", "orange", "lightpink")
plot(predicted_lulc_2035, main = "Predicted LULC 2035", col = colors)
legend("topright", legend = lulc_class_names, fill = colors, bty = "n")




# VISUALIZATION


# --- 1. Prepare Data Frames for Each Predicted Year ---

# Define your LULC classes (ensure this matches your previous setup)
lulc_class_names <- c("Water", "Forest", "Built-up", "Aquaculture", "Others")
num_lulc_classes <- length(lulc_class_names)

# Define a consistent color palette for all plots
custom_colors <- c("Water" = "lightblue",
                   "Forest" = "darkgreen",
                   "Built-up" = "orange",
                   "Aquaculture" = "darkblue",
                   "Others" = "red")

# Function to prepare a single predicted raster for ggplot
prepare_for_ggplot <- function(raster_obj, year_label) {
  df <- as.data.frame(raster_obj, xy = TRUE, na.rm = TRUE)
  colnames(df)[3] <- "LULC_Class"
  df$LULC_Class <- factor(df$LULC_Class,
                          levels = 0:(num_lulc_classes - 1),
                          labels = lulc_class_names)
  df$Year <- factor(year_label) # Add a year column for consistency if needed later, though not strictly for side-by-side
  return(df)
}

# Convert each predicted raster to a ggplot-ready data frame
lulc_2021 <- raster("Data/sundarban tifs/Sundarban_rf_21.tif")
lulc_df_2021 <- prepare_for_ggplot(lulc_2021, "2021")
lulc_df_2025 <- prepare_for_ggplot(predicted_lulc_2025, "2025")
lulc_df_2030 <- prepare_for_ggplot(predicted_lulc_2030, "2030")
lulc_df_2035 <- prepare_for_ggplot(predicted_lulc_2035, "2035")

# --- 2. Create Individual ggplot Objects ---

# Plot for 2025
p_2021 <- ggplot() +
  geom_raster(data = lulc_df_2021, aes(x = x, y = y, fill = LULC_Class)) +
  scale_fill_manual(values = custom_colors, name = "LULC Class") +
  labs(title = "Observed LULC: 2021", 
       x = "Longitude", 
       y = "Latitude",
       fill = "2021 LULC Class") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), # Center title
        axis.text = element_blank(), # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        axis.title = element_blank()) + # Remove axis titles
  coord_sf(expand = FALSE)

p_2025 <- ggplot() +
  geom_raster(data = lulc_df_2025, aes(x = x, y = y, fill = LULC_Class)) +
  scale_fill_manual(values = custom_colors, name = "LULC Class") +
  labs(title = "Predicted LULC: 2025", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), # Center title
        axis.text = element_blank(), # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        axis.title = element_blank()) + # Remove axis titles
  coord_sf(expand = FALSE) # Important for spatial plots, remove extra whitespace

# Plot for 2030
p_2030 <- ggplot() +
  geom_raster(data = lulc_df_2030, aes(x = x, y = y, fill = LULC_Class)) +
  scale_fill_manual(values = custom_colors, name = "LULC Class") +
  labs(title = "Predicted LULC: 2030", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") + # Hide legend for middle plot
  coord_sf(expand = FALSE)

# Plot for 2035
p_2035 <- ggplot() +
  geom_raster(data = lulc_df_2035, aes(x = x, y = y, fill = LULC_Class)) +
  scale_fill_manual(values = custom_colors, name = "LULC Class") +
  labs(title = "Predicted LULC: 2035", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  coord_sf(expand = FALSE)


# --- 3. Arrange Plots Side-by-Side using patchwork ---

# Use '+' to add plots side-by-side
# You can also use '/' for stacking, or specify layout with plot_layout()
combined_plot <- p_2021 + p_2030  +
  plot_layout(guides = "collect")  # Position the collected legend

# Display the combined plot
print(combined_plot)

# --- Optional: Save the combined plot ---
ggsave("predicted_lulc_2025_2030_2035_side_by_side.png", 
       plot = combined_plot, 
       width = 15, 
       height = 6, 
       dpi = 300
       )

