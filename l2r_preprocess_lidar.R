# Preprocess lidar for L2R (Webster et al., 2020)
# Nic Tarasewicz
# 07-10-2025

# Run for each area. Un-comment for checks. 
# Update i/o PATHS before running 

# Load libraries
library(lidR)
library(sf)
library(dplyr)
library(data.table)
library(terra)
library(rlas)

# Explicitly use lidR's crs functions to avoid masking issues
lidR_crs <- lidR::crs
lidR_crs_set <- lidR::`crs<-`

# NEON lidar file paths
file_paths <- c(
  "*/NEON_provisional_2023/NEON.D13.NIWO.DP1.30003.001.2023-08.basic.20241202T020119Z.PROVISIONAL/NEON_D13_NIWO_DP1_452000_4431000_classified_point_cloud_colorized.laz",
  "*/NEON_provisional_2023/NEON.D13.NIWO.DP1.30003.001.2023-08.basic.20241202T020119Z.PROVISIONAL/NEON_D13_NIWO_DP1_452000_4432000_classified_point_cloud_colorized.laz",
  "*/NEON.D13.NIWO.DP1.30003.001.2023-08.basic.20241202T020119Z.PROVISIONAL/NEON_D13_NIWO_DP1_453000_4431000_classified_point_cloud_colorized.laz",
  "*/NEON_provisional_2023/NEON.D13.NIWO.DP1.30003.001.2023-08.basic.20241202T020119Z.PROVISIONAL/NEON_D13_NIWO_DP1_453000_4432000_classified_point_cloud_colorized.laz"
)

# CRS
desired_crs <- "+proj=utm +zone=13 +datum=WGS84 +units=m +no_defs"

# Load .laz
las_list <- lapply(file_paths, function(file) {
  las <- readLAS(file)
  if (is.null(lidR_crs(las)) || is.na(lidR_crs(las))) {
    lidR_crs_set(las) <- desired_crs  # Manually set the CRS
  }
  return(las)
})
merged_las <- do.call(rbind, las_list)

# Clip to domain area
# area_shape <- st_read("*/forest_fen_area.shp") # EcoTram domain 
# area_shape <- st_read("*/continuous_forest.shp") # US-NR1 domain
clipped_las <- clip_roi(merged_las, area_shape)

### Cleaning provisional NEON data
las_check(clipped_las)

# Remove duplicated points
clipped_las <- filter_duplicates(clipped_las)

# Handle points with identical X, Y but differing Z
degenerated_points <- filter_poi(clipped_las, Classification == LASGROUND)
unique_xy <- unique(degenerated_points@data[, c("X", "Y")])
clipped_las <- filter_poi(clipped_las, !duplicated(paste(X, Y)))

# Run las_check again to confirm all issues are resolved
las_check(clipped_las)

# Classify and remove noise
clipped_las <- classify_noise(clipped_las, ivf(5,2))
clipped_las <- filter_poi(clipped_las, Classification != LASNOISE)

# Classify ground points
clipped_las <- classify_ground(clipped_las, csf(sloop_smooth = TRUE, 
                                                class_threshold = 0.5, 
                                                cloth_resolution = 1, 
                                                time_step = 1))

# Remove points classified as 1 (Unclassified) and 7 (Noise)
clipped_las <- filter_poi(clipped_las, Classification != 1)
clipped_las <- filter_poi(clipped_las, Classification != 7)
table(clipped_las$Classification)
plot(clipped_las, color = "Classification")
writeLAS(clipped_las, "DOMAIN lidar.laz")

### Create DTM
dtm <- grid_terrain(clipped_las, res = 0.5, algorithm = tin())
Sys.setenv("GDAL_PAM_ENABLED" = "YES")
writeRaster(dtm, "Domain DTM PATH.tif", overwrite = TRUE,
            wopt = list(gdal = c("TFW=YES")))


### Create CHM (pit-free)
normalized_las <- normalize_height(clipped_las, tin())
plot(normalized_las, color = "Classification")
writeLAS(normalized_las, "Domain CHM PATH.laz") # Necessary for L2R input

chm <- rasterize_canopy(
  normalized_las,
  res = 0.5,
  pitfree(
    thresholds = c(0, 10, 20),  # Sequential height thresholds
    max_edge = c(0, 1.5, 1.5)   # Max edge lengths for TIN
  )
)
plot(chm, col = terrain.colors(25))
writeRaster(chm, "Domain CHM PATH.tif", overwrite = TRUE,
            wopt = list(gdal = c("TFW=YES")))


### Create DSM (pit-free)
dsm_pitfree <- rasterize_canopy(
  clipped_las,
  res = 0.5,
  pitfree(thresholds = c(0, 10, 20), max_edge = c(0, 1.5, 1.5))
)
plot(dsm_pitfree, col = terrain.colors(25))
writeRaster(dsm_pitfree, "Domain DSM PATH.tif", overwrite = TRUE,
            wopt = list(gdal = c("TFW=YES")))


### Create DTM buffer for wider domain area (necessary for L2R input)
area_shape_buffer <- st_read("*/niwot_buffer.shp")
clipped_las_buffer <- clip_roi(merged_las, area_shape_buffer)

# Remove duplicated points
clipped_las_buffer <- filter_duplicates(clipped_las_buffer)
# Handle points with identical X, Y but differing Z
degenerated_points <- filter_poi(clipped_las_buffer, Classification == LASGROUND)
unique_xy <- unique(degenerated_points@data[, c("X", "Y")])
clipped_las_buffer <- filter_poi(clipped_las_buffer, !duplicated(paste(X, Y)))

# Classify and remove noise
clipped_las_buffer <- classify_noise(clipped_las_buffer, ivf(5,2))
clipped_las_buffer <- filter_poi(clipped_las_buffer, Classification != LASNOISE)

# Classify ground points
# Ensure ground points are classified (if not already)
clipped_las_buffer <- classify_ground(clipped_las_buffer, csf(sloop_smooth = TRUE, 
                                                class_threshold = 0.5, 
                                                cloth_resolution = 1, 
                                                time_step = 1))

# Remove points classified as 1 (Unclassified) and 7 (Noise)
clipped_las_buffer <- filter_poi(clipped_las_buffer, Classification != 1)
clipped_las_buffer <- filter_poi(clipped_las_buffer, Classification != 7)
table(clipped_las_buffer$Classification)
plot(clipped_las_buffer, color = "Classification")

# Create DTM
dtm_buff <- grid_terrain(clipped_las_buffer, res = 0.5, algorithm = tin())
plot(dtm_buff)
writeRaster(dtm_buff, "DTM buffer PATH.tif", overwrite = TRUE,
            wopt = list(gdal = c("TFW=YES")))


### Tree segmentation
las <- readLAS("NORMALIZED LIDAR .LAZ PATH", filter = "-drop_z_below 0")
chm <- rasterize_canopy(las, 0.5, pitfree())

ttops <- locate_trees(las, lmf(ws = 3, hmin = 2)) # identify tree tops

# Convert tree tops to a 2D spatial object
ttops_sf <- st_as_sf(ttops, coords = c("x", "y"), crs = st_crs(chm))
ttops_2D <- st_zm(ttops_sf, drop = TRUE, what = "ZM")

st_write(ttops_2D, "tree tops PATH.shp")

# Visual check
plot(chm, col = height.colors(50))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)
x <- plot(las, bg = "white", size = 4)
add_treetops3d(x, ttops)

# Segmentation of point-cloud
algo <- dalponte2016(chm, ttops)
crowns <- segment_trees(las, algo)
plot(crowns, bg = "white", size = 4, color = "treeID")

# Filter point cloud
filtered_las <- filter_poi(crowns, !is.na(treeID)) 
plot(filtered_las, bg = "white", size = 4, color = "treeID")
plot(filtered_las, color = "Classification")
table(filtered_las$Classification)
writeLAS(filtered_las, "tree filtered PATH.laz")

# Adjust lidar to match (required by L2R)
filtered_las@data$Intensity <- as.integer(filtered_las@data$treeID)
# Rescale and reoffset
xmin <- min(filtered_las@data$X)
ymin <- min(filtered_las@data$Y)
zmin <- min(filtered_las@data$Z)
filtered_las <- las_rescale(filtered_las, 0.001, 0.001, 0.001)
filtered_las <- las_reoffset(filtered_las, xmin, ymin, zmin)
writeLAS(filtered_las, "FILTERED POINT CLOUD PATH.laz") # necessary for L2R input


### Tree metrics
filtered_las <- readLAS("FILTERED POINT CLOUD PATH.laz")
# Crown polygons
crown_polygons <- crown_metrics(filtered_las, func = .stdtreemetrics, geom = "convex")
st_write(crown_polygons, "TREE CROWNS PATH.shp")
str(crown_polygons)
st_bbox(crown_polygons)
colnames(crown_polygons)
summary(filtered_las)

# Compute new tree tops as centroids of crown polygons
ttops <- st_centroid(crown_polygons)
st_writest_write(
  ttops,
  "TREE TOP CENTROIDS.shp",
  delete_layer = TRUE
)

# Compute DBH for each tree (Jucker et al., 2017 - Biome 10: Nearctic Temperate coniferous Gymnosperm)
alpha <- 0.589
beta  <- 0.817

cd_m <- 2 * sqrt(crown_polygons$convhull_area / pi)  # diameter in meters
crown_polygons$DBH <- alpha * (crown_polygons$Z * cd_m)^beta  # DBH in cm

# Combine centroids (tree tops) with metrics
ttops_with_metrics <- cbind(ttops, crown_polygons[, c("treeID", "Z", "DBH")])

# Generate .txt file
output_txt <- data.table(
  treeptsX = st_coordinates(ttops_with_metrics)[, 1],  # X coordinate of centroid
  treeptsY = st_coordinates(ttops_with_metrics)[, 2],  # Y coordinate of centroid
  treeptsH_m = ttops_with_metrics$Z,                  # Tree height (Z)
  treeptsdbh_cm = ttops_with_metrics$DBH,             # DBH
  tree_num = ttops_with_metrics$treeID                # Tree ID
)
# Sanity filters
n_before <- nrow(output_txt)
output_txt <- output_txt[!is.na(treeptsH_m) & !is.na(treeptsdbh_cm) & treeptsH_m > 0 & treeptsdbh_cm > 0]
n_after <- nrow(output_txt)
if (n_after < n_before) {
  message(sprintf("Filtered %d rows with invalid height/DBH.", n_before - n_after))
}
fwrite(output_txt, "treeinfo PATH.txt", sep = "\t") # necessary for L2R input


### Generate points from .shp for L2R model locations within domain
process_shp_to_txt <- function(shp_file_path, output_txt_file_path) {
  shp_data <- st_read(shp_file_path, quiet = TRUE)
  if (!all(st_geometry_type(shp_data) == "POINT")) stop("The .shp file must contain only POINT geometries.")
  coords <- st_coordinates(shp_data)
  fwrite(
    data.table(X = coords[, 1], Y = coords[, 2]),
    output_txt_file_path, sep = "\t", col.names = FALSE, row.names = FALSE
  )
  cat("Output .txt file generated at:", output_txt_file_path, "\n")
}

shp_file_path <- "DOMAIN AREA POINTS.shp"
output_txt_file_path <- "DOMAIN POINTS.txt"
process_shp_to_txt(shp_file_path, output_txt_file_path)