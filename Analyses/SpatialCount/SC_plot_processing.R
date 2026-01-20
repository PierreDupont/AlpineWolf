## -------------------------------------------------------------------------- ##
## ------------------------ ALPINE WOLF SC ---------------------------------- ##
## -------------------------------------------------------------------------- ##
## ------ CLEAN THE WORK ENVIRONMENT ------

rm(list=ls())

## ------ IMPORT REQUIRED LIBRARIES ------
library(raster)
library(coda)
library(nimble)
library(nimbleSCR)
library(stringr)
library(abind)
library(R.utils)
library(sf)
library(fasterize)
library(dplyr)
library(tidyr)
library(lubridate)
library(stars)
library(units)
library(RANN)
library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
library(gridExtra)
library(MetBrewer)
library(data.table)
library(ggplot2)
library(fields)
library(cowplot)    # for arranging multiple plots



## ------ SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)
sourceCpp(file = file.path(getwd(),"Source/cpp/GetDensity.cpp"))
sourceCpp(file = file.path(getwd(),"Source/cpp/GetSpaceUse.cpp"))

## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
modelName = "AlpineWolf.SC.Plots"
thisDir <- file.path(analysisDir, modelName)
# setwd("/Users/virginia/Dropbox/AlpineWolf/02_Analysis/AlpineWolf.SC.Plots")

if(is.null(modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
if(!dir.exists(thisDir)){dir.create(thisDir)}

## HABITAT SPECIFICATIONS
habitat = list( resolution = 5000, 
                buffer = 30000,
                dOK = 100, 
                dmax = 25)

## ------ I. LOAD & CLEAN DATA ------
## ------   1. STUDY AREA ------
##---- Polygon of Italy and neighbouring countries
countries <- read_sf(file.path(dataDir,"GISData/Countries/Europe_32N.shp"))
countries <- st_transform(countries, 32632)

##--- Study area grid
studyAreaGrid <- read_sf(file.path(dataDir,"GISData/shape_studyarea_ALPS/Studyarea_ALPS_2020_2021.shp"))
studyAreaGrid <- st_transform(x = studyAreaGrid, crs = st_crs(countries))
studyArea <- studyAreaGrid %>%
  st_snap(x = ., y = ., tolerance = 0.0001) %>%
  st_union() 

##---- create raster with desired resolution
grid.r <- st_as_stars(st_bbox(studyArea),
                      dx = habitat$resolution,
                      dy = habitat$resolution)
##---- mask and crop
grid.r <- grid.r[studyArea, ] 
##---- convert from stars to sf objects
grid <- st_as_sf(grid.r)
grid <- grid[studyArea, ] 
##---- give an ID to each grid cell 
grid$id <- 1:nrow(grid)
##---- Remove the "values" column
grid <- grid[-1]

##---- create centroids 
grid$centroids <- st_centroid(grid) %>% 
  st_geometry() 

## ------   2. HABITAT 
## ------     2.1. HABITAT CHARACTERISTICS ------
##---- Create Habitat polygon (detector grid + buffer)
habitat$polygon <- st_buffer(st_union(grid),
                             habitat$buffer)

##---- Remove un-suitable habitat (e.g. seas)
habitat$polygon <- st_intersection( st_union(countries),
                                    habitat$polygon,
                                    drop_lower_td = TRUE)

##---- Create habitat raster (keep raster cells with > 50% habitat)
##---- (the huge buffer is used to make sure all models have comparable rasters)
habitat$raster <- raster(extent(st_bbox(st_buffer(st_union(studyAreaGrid),
                                                  100000))))

##---- Trick to get % of habitat
res(habitat$raster) <- habitat$resolution/10
habitat$raster <- fasterize( countries,
                             habitat$raster, )
habitat$raster[is.na(habitat$raster)]<- 0
habitat$raster <- aggregate(habitat$raster, fact = 10, sum, na.rm = TRUE)
habitat$raster[habitat$raster[ ] < 50] <- NA
habitat$raster[habitat$raster[ ] >= 50] <- 1
# plot(habitat$raster)


##---- Mask and crop habitat to the buffer
habitat$raster <- mask(habitat$raster, st_as_sf(habitat$polygon))
habitat$raster <- crop(habitat$raster, st_as_sf(habitat$polygon))

##---- Create habitat grid
habitat$grid <- st_as_sf(rasterToPolygons( habitat$raster,
                                           fun = function(x)x == 1))
habitat$grid$id <- 1:nrow(habitat$grid)
habitat$grid <- habitat$grid[,c(2,3)]
st_crs(habitat$grid) <- st_crs(habitat$polygon)

##---- Create Habitat matrix of cell ID
habitat$matrix <- habitat$raster
habitat$matrix[] <- 0

##----  Identify suitable habitat cells
isHab <- which(habitat$raster[]==1)

##---- Cell ID starts from the top left corner, increments left to right and
##---- up to down
habitat$matrix[isHab] <- 1:length(isHab)

##---- Convert to matrix
habitat$matrix <- as.matrix(habitat$matrix)
habitat$binary <- habitat$matrix  ##(required by the getLocalObject function; the function could be changed to deal with matrix others than 0/1)
habitat$binary[habitat$binary > 0] <- 1

##---- Obtain xy coordinates of habitat cells
habitat$coords <- coordinates(habitat$raster)[isHab, ]
dimnames(habitat$coords) <- list(1:length(isHab), c("x","y"))
habitat$sp <- SpatialPoints(coords = habitat$coords,
                            proj4string = crs(habitat$polygon))

##---- Retrieve habitat windows corners
habitat$lowerCoords <- habitat$coords - 0.5*habitat$resolution
habitat$upperCoords <- habitat$coords + 0.5*habitat$resolution
habitat$n.HabWindows <- dim(habitat$lowerCoords)[1] ## == length(isHab)


##----  Extract the habitat raster
habitat.r <- habitat$raster

##---- Set cells outside the habitat to NA
habitat.r[habitat.r == 0] <- NA

##---- Create a matrix of raster cellsIDs
habitat.id <- matrix( data = 1:ncell(habitat.r),
                      nrow = dim(habitat.r)[1],
                      ncol = dim(habitat.r)[2],
                      byrow = TRUE)
##---- Create "ITALY" habitat raster (for density estimates)
habitat$Italia <- mask(habitat$raster, st_as_sf(studyArea))
habitat$Italia[habitat$Italia == 0] <- NA

## ------   3. LOAD MODELS POSTERIORS ------

## Step 1: List only the 5 model files (exclude the reference one)
# Use pattern matching to exclude the reference model file
modelFiles <- list.files(thisDir,pattern = "_posterior\\.RData$")
modelFiles <- modelFiles[!grepl("AlpineWolf.5.2_RJMCMC", modelFiles)]  # adjust as needed

## Step 2: Load WA_Italy from the 5 new models
model_outputs <- list()
for (i in seq_along(modelFiles)) {
  fullPath <- file.path(thisDir, modelFiles[i])  # construct full path
  load(fullPath)                                 # load the RData file
  model_outputs[[i]] <- WA_Italy                 # store WA_Italy into list
}

## Step 3: Load reference model (WA_regions)
ref_file <- "AlpineWolf.5.2_RJMCMC_posterior.RData"  # your reference filename
load(file.path(thisDir, ref_file))                   # correct usage
model_outputs[[6]] <- WA_regions

## Step 4: Set names
names(model_outputs) <- c(paste0("M", 1:5), "SCR Marucco et al. 2023")

# Get minimum chain length
min_iter <- min(sapply(model_outputs, function(x) length(x$PosteriorAllRegions)))

# Truncate all to same length
df <- data.frame(
  M1  = model_outputs[[1]][["PosteriorAllRegions"]][1:min_iter],
  M2  = model_outputs[[2]][["PosteriorAllRegions"]][1:min_iter],
  M3  = model_outputs[[3]][["PosteriorAllRegions"]][1:min_iter],
  M4  = model_outputs[[4]][["PosteriorAllRegions"]][1:min_iter],
  M5  = model_outputs[[5]][["PosteriorAllRegions"]][1:min_iter],
  SCR = model_outputs[[6]][["PosteriorAllRegions"]][1:min_iter]
)
  
# Convert to long format
df_long <- df %>%
  tidyr::pivot_longer(cols = everything(), names_to = "model", values_to = "N")

df_long <- data.frame(
  N = model_outputs[[1]][["PosteriorAllRegions"]],
  model = "M1"
)
for (i in 2:6) {
  df_long <- rbind(
    df_long,
    data.frame(
      N = model_outputs[[i]][["PosteriorAllRegions"]],
      model = c("M2", "M3", "M4", "M5", "SCR")[i - 1]
    )
  )
}
## ------ 4 DENSITY MAPS FOR ALL MODELS -------
## ------     4.1 R PLOT VERSION -----------
# tiff("/Users/virginia/Desktop/WolfDensity_2x3_LegendRight.tiff", 
#      width = 2400, height = 1600, res = 300)
#   
# # This creates a mask of where data is not NA
# nonNA_mask <- habitat$Italia
# nonNA_mask[] <- ifelse(!is.na(habitat$Italia[]), 1, NA)
# extent_mask <- extent(nonNA_mask)
# 
# # Add a slight buffer
# buffer <- 20000
# xlim_use <- c(xmin(extent_mask) - buffer, xmax(extent_mask) + buffer)
# ylim_use <- c(ymin(extent_mask) - buffer, ymax(extent_mask) + buffer)
# 
# # Crop countries to match
# countries_crop <- st_crop(countries, xmin = xlim_use[1], xmax = xlim_use[2],
#                           ymin = ylim_use[1], ymax = ylim_use[2])
# 
# # Step 1: Convert country polygons to boundaries
# boundaries <- st_boundary(countries_crop)
# 
# # Step 2: Wrap as geometry column and cast to MULTILINESTRING
# boundary_lines <- st_cast(st_sf(geometry = boundaries), "MULTILINESTRING")
# 
# # Step 3: Done — no need to filter manually
# countries_lines <- boundary_lines
# # Color and page setup
# layout_matrix <- matrix(c(1, 2, 3, 7,
#                           4, 5, 6, 7), nrow = 2, byrow = TRUE)
# 
# layout(mat = layout_matrix, widths = c(1, 1, 1, 0.5), heights = c(1, 1))
# par(mfrow = c(2, 3), mar = c(2, 2, 2, 2))
# maxDens <- max(sapply(model_outputs, function(x) max(x$MeanCell, na.rm = TRUE)))
# cuts <- seq(0, maxDens, length.out = 100)
# colFunc <- colorRampPalette(c("white", "skyblue", "steelblue3", "lightgreen", "coral", "yellow"))
# col <- colFunc(length(cuts) - 1)
# 
# for (i in 1:6) {
#   dens.R <- raster(habitat$Italia)
#   dens.R[] <- model_outputs[[i]]$MeanCell
#   dens.R[is.na(habitat$Italia[])] <- NA
#   
#   plot(dens.R,
#        breaks = cuts, col = col,
#        legend = FALSE,
#        xlim = xlim_use, ylim = ylim_use,
#        axes = FALSE, box = FALSE)
#   
#   plot(habitat$polygon, add = TRUE, border = grey(0.6), lwd = 1)
#   plot(st_geometry(countries_lines), add = TRUE, border = "black", lwd = 0.5)
#   
#   title(main = names(model_outputs)[i], font.main = 2, cex.main = 1.2, line = 0.2)
#   N_est <- round(model_outputs[[i]]$summary["Total", 1], 1)
#   N_lwr <- round(model_outputs[[i]]$summary["Total", 4], 1)
#   N_upr <- round(model_outputs[[i]]$summary["Total", 5], 1)
#   mtext(paste0("N = ", N_est, " [", N_lwr, " ; ", N_upr, "]"),
#         side = 1, line = 0.5, cex = 0.9)
# }
# 
# # Now manually place the legend using par(fig=...) — full control
# par(mar = c(4, 1, 4, 4))  # wider right margin
# image.plot(legend.only = TRUE,
#            zlim = c(0, maxDens),
#            col = col, breaks = cuts,
#            legend.width = 2,
#            legend.shrink = 0.8,
#            axis.args = list(cex.axis = 1.2),
#            legend.args = list(text = "Wolf Density", side = 3, line = 1, cex = 1.2))
# 
# 
# 
# dev.off()
# 
# 
# 
# # Save the output image
# tiff("/Users/virginia/Desktop/WolfDensity_2x3_LegendRight.tiff", 
#      width = 2400, height = 1600, res = 300)
# 
# # Create a 2x4 layout: 6 maps and 1 legend panel
# layout_matrix <- matrix(c(1, 2, 3, 7,
#                           4, 5, 6, 7), nrow = 2, byrow = TRUE)
# layout(mat = layout_matrix, widths = c(1, 1, 1, 0.5), heights = c(1, 1))
# 
# # Adjust margins
# par(mar = c(1.5, 1.5, 2, 1))
# 
# # Get common scale
# maxDens <- max(sapply(model_outputs, function(x) max(x$MeanCell, na.rm = TRUE)))
# cuts <- seq(0, maxDens, length.out = 100)
# colFunc <- colorRampPalette(c("white", "skyblue", "steelblue3", "lightgreen", "coral", "yellow"))
# col <- colFunc(length(cuts) - 1)
# 
# # Plot each model
# for (i in 1:6) {
#   dens.R <- raster(habitat$Italia)
#   dens.R[] <- model_outputs[[i]]$MeanCell
#   dens.R[is.na(habitat$Italia[])] <- NA
#   
#   plot(dens.R,
#        breaks = cuts, col = col,
#        legend = FALSE,
#        xlim = xlim_use, ylim = ylim_use,
#        axes = FALSE, box = FALSE)
#   
#   plot(habitat$polygon, add = TRUE, border = grey(0.6), lwd = 1)
#   plot(st_geometry(countries_lines), add = TRUE, border = "black", lwd = 0.5)
#   
#   title(main = names(model_outputs)[i], font.main = 2, cex.main = 1.1, line = 0.5)
#   N_est <- round(model_outputs[[i]]$summary["Total", 1], 1)
#   N_lwr <- round(model_outputs[[i]]$summary["Total", 4], 1)
#   N_upr <- round(model_outputs[[i]]$summary["Total", 5], 1)
#   mtext(paste0("N = ", N_est, " [", N_lwr, " ; ", N_upr, "]"),
#         side = 1, line = 0.5, cex = 0.8)
# }
# 
# # Legend panel
# par(mar = c(3, 1, 3, 3))
# image.plot(legend.only = TRUE,
#            zlim = c(0, maxDens),
#            col = col, breaks = cuts,
#            legend.width = 1.5,
#            legend.shrink = 0.7,
#            axis.args = list(cex.axis = 1),
#            legend.args = list(text = "Wolf Density", side = 3, line = 0.5, cex = 1))
# 
# dev.off()

## ------     4.2 GGPLOT VERSION -----------  

# Convert vector data to sf (if not already) and ensure they share the raster CRS
habitat_sf    <- st_as_sf(habitat$polygon)               # study area outline as sf
countries_sf  <- st_as_sf(countries)               # country borders as sf
# Transform country borders to the habitat's coordinate system if needed
countries_sf  <- st_transform(countries_sf, st_crs(habitat_sf))

# Determine a common color scale range across all six density rasters
# Create list of density rasters (masked by habitat)
density_rasters <- lapply(model_outputs, function(x) {
  r <- raster(habitat$Italia)
  crs(r) <- st_crs(habitat_sf)$proj4string
  r[] <- x$MeanCell
  r[is.na(habitat$Italia[])] <- NA
  return(r)
})

# Assign to named objects for convenience (optional)
density_raster1 <- density_rasters[[1]]
density_raster2 <- density_rasters[[2]]
density_raster3 <- density_rasters[[3]]
density_raster4 <- density_rasters[[4]]
density_raster5 <- density_rasters[[5]]
density_raster6 <- density_rasters[[6]]


# Now calculate global min and max for color scale
global_min <- min(sapply(density_rasters, function(r) raster::cellStats(r, stat = 'min')), na.rm = TRUE)
global_max <- max(sapply(density_rasters, function(r) raster::cellStats(r, stat = 'max')), na.rm = TRUE)

# Convert each raster to a data frame of pixel coordinates and values for ggplot2
df1 <- as.data.frame(density_raster1, xy = TRUE)
names(df1) <- c("x", "y", "density")
df1 <- df1[!is.na(df1$density), ]  # optional: remove NA cells for efficiency

df2 <- as.data.frame(density_raster2, xy = TRUE)
names(df2) <- c("x", "y", "density")
df2 <- df2[!is.na(df2$density), ]

df3 <- as.data.frame(density_raster3, xy = TRUE)
names(df3) <- c("x", "y", "density")
df3 <- df3[!is.na(df3$density), ]

df4 <- as.data.frame(density_raster4, xy = TRUE)
names(df4) <- c("x", "y", "density")
df4 <- df4[!is.na(df4$density), ]

df5 <- as.data.frame(density_raster5, xy = TRUE)
names(df5) <- c("x", "y", "density")
df5 <- df5[!is.na(df5$density), ]

df6 <- as.data.frame(density_raster6, xy = TRUE)
names(df6) <- c("x", "y", "density")
df6 <- df6[!is.na(df6$density), ]

# Define model names and corresponding abundance (N) estimates with CIs for captions
model_names <- c("M1", "M2", "M3", "M4", "M5", "SCR Marucco et al. 2023")
# Extract abundance point estimates and 95% CIs from model summaries
N_est   <- sapply(model_outputs, function(x) round(x$summary["Total", 1], 1))  # mean
N_lower <- sapply(model_outputs, function(x) round(x$summary["Total", 4], 1))  # 2.5%
N_upper <- sapply(model_outputs, function(x) round(x$summary["Total", 5], 1))  # 97.5%

# Determine common spatial extent (bounding box) to use for all plots (ensures same zoom)
plot_bbox <- st_bbox(studyArea)  # use study area extent
xlim <- c(plot_bbox["xmin"], plot_bbox["xmax"])
ylim <- c(plot_bbox["ymin"], plot_bbox["ymax"])

# Define custom color palette and breaks
custom_colors <- c("white","steelblue3", "skyblue", "#6fd1be", "lightgreen", "yellow")
n_breaks <- 100
col_palette <- colorRampPalette(custom_colors)(n_breaks)

# Create ggplot maps for each model
p1 <- ggplot() +
  geom_raster(data = df1, aes(x = x, y = y, fill = density)) +
  geom_sf(data = countries_sf, fill = NA, color = "black", size = 0.1) +
  # geom_sf(data = habitat_sf, fill = NA, color = "gray40", size = 1) +
  scale_fill_gradientn(colors = col_palette, 
                       limits = c(global_min, global_max),
                       name = "Wolf density",
                       na.value = NA) +
  labs(title = model_names[1],
       caption = paste0("N = ", N_est[1], " [", N_lower[1], " ; ", N_upper[1], "]")) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_void() +  # removes all axes, ticks, and background grid
  theme(
    legend.position = "none",
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5)
  )

p2 <- ggplot() +
  geom_raster(data = df2, aes(x = x, y = y, fill = density)) +
  geom_sf(data = countries_sf, fill = NA, color = "black", size = 0.1) +
  # geom_sf(data = habitat_sf, fill = NA, color = "gray40", size = 1) +
  scale_fill_gradientn(colors = col_palette, 
                       limits = c(global_min, global_max),
                       name = "Wolf density",
                       na.value = NA) +
  labs(title = model_names[2],
       caption = paste0("N = ", N_est[2], " [", N_lower[2], " ; ", N_upper[2], "]")) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_void() +  # removes all axes, ticks, and background grid
  theme(
    legend.position = "none",
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5)
  )

p3 <- ggplot() +
  geom_raster(data = df3, aes(x = x, y = y, fill = density)) +
  geom_sf(data = countries_sf, fill = NA, color = "black", size = 0.1) +
  # geom_sf(data = habitat_sf, fill = NA, color = "gray40", size = 1) +
  scale_fill_gradientn(colors = col_palette, 
                       limits = c(global_min, global_max),
                       name = "Wolf density",
                       na.value = NA) +
  labs(title = model_names[3],
       caption = paste0("N = ", N_est[3], " [", N_lower[3], " ; ", N_upper[3], "]")) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_void() +  # removes all axes, ticks, and background grid
  theme(
    legend.position = "none",
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5)
  )

p4 <- ggplot() +
  geom_raster(data = df4, aes(x = x, y = y, fill = density)) +
  geom_sf(data = countries_sf, fill = NA, color = "black", size = 0.1) +
  # geom_sf(data = habitat_sf, fill = NA, color = "gray40", size = 1) +
  scale_fill_gradientn(colors = col_palette, 
                       limits = c(global_min, global_max),
                       name = "Wolf density",
                       na.value = NA) +
  labs(title = model_names[4],
       caption = paste0("N = ", N_est[4], " [", N_lower[4], " ; ", N_upper[4], "]")) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_void() +  # removes all axes, ticks, and background grid
  theme(
    legend.position = "none",
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5)
  )

p5 <- ggplot() +
  geom_raster(data = df5, aes(x = x, y = y, fill = density)) +
  geom_sf(data = countries_sf, fill = NA, color = "black", size = 0.1) +
  # geom_sf(data = habitat_sf, fill = NA, color = "gray40", size = 1) +
  scale_fill_gradientn(colors = col_palette, 
                       limits = c(global_min, global_max),
                       name = "Wolf density",
                       na.value = NA) +
  labs(title = model_names[5],
       caption = paste0("N = ", N_est[5], " [", N_lower[5], " ; ", N_upper[5], "]")) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_void() +  # removes all axes, ticks, and background grid
  theme(
    legend.position = "none",
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5)
  )

p6 <- ggplot() +
  geom_raster(data = df6, aes(x = x, y = y, fill = density)) +
  geom_sf(data = countries_sf, fill = NA, color = "black", size = 0.1) +
  # geom_sf(data = habitat_sf, fill = NA, color = "gray40", size = 1) +
  scale_fill_gradientn(colors = col_palette, 
                       limits = c(global_min, global_max),
                       name = "Wolf density",
                       na.value = NA) +
  labs(title = model_names[6],
       caption = paste0("N = ", N_est[6], " [", N_lower[6], " ; ", N_upper[6], "]")) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_void() +  # removes all axes, ticks, and background grid
  theme(
    legend.position = "none",
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5)
  )

# Use cowplot to arrange the six maps and add a shared legend on the right
# Extract a legend from one plot (use the first plot with legend on to get the fill scale legend)
# legend_plot <- ggplot(df1, aes(x = x, y = y, fill = density)) +
#   geom_raster() +
#   scale_fill_gradientn(
#     colors = col_palette,
#     limits = c(global_min, global_max),
#     name = "Wolf density",
#     na.value = NA
#   ) +
#   theme_minimal() +
#   theme(
#     legend.position = "bottom",
#     legend.direction = "horizontal",
#     legend.title = element_text(size = 11),
#     legend.text = element_text(size = 10),
#     legend.key.width = unit(2, "cm"),
#     legend.key.height = unit(0.4, "cm"),
#     plot.margin = margin(0, 0, 0, 0)
#   )


legend_plot <- ggplot() + geom_raster(data = df1, aes(x = x, y = y, fill = density)) + 
  scale_fill_gradientn(colors = col_palette, 
                      limits = c(global_min, global_max), 
                      name = "Wolf density", na.value = NA) + 
  theme_minimal() + 
  theme(legend.position = "right", 
        legend.title = element_text(size = 11), 
        legend.text = element_text(size = 10), 
        plot.background = element_blank()) # no plot, just legend 
legend_grob <- cowplot::get_legend(legend_plot)

# Arrange plots in a 2x3 grid (no legends in these plots)
plot_grid_2x3 <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, 
                                    legend_grob,
                                    ncol = 3, nrow = 2)

# Combine the grid and legend (now below)
combined_plot <- cowplot::plot_grid(
  plot_grid_2x3,
  legend_grob, 
  ncol = 2, 
  rel_widths = c(1, 0.12))
# Save the combined figure to a TIFF file at 2400x1600 px, 300 dpi
ggsave("/Users/virginia/Desktop/WolfDensity_2x3_LegendRight.png",
       plot = combined_plot, 
       width = 2400/300, 
       height = 1600/300, 
       dpi = 300,
       units = "in",
       device = "png")

# ggsave("/Users/virginia/Desktop/WolfDensity_2x3_LegendBottom.png", 
#        plot = combined_plot, 
#        width = 2400 / 300, 
#        height = 1700 / 300,  # slightly taller to fit legend
#        dpi = 300, 
#        units = "in",
#        device = "png")

## ------   5. VIOLIN PLOTS OF POSTERIORS ------

N_all <- rbind(
  data.frame(N = model_outputs[[1]][["PosteriorAllRegions"]], model = "M1"),
  data.frame(N = model_outputs[[2]][["PosteriorAllRegions"]], model = "M2"),
  data.frame(N = model_outputs[[3]][["PosteriorAllRegions"]], model = "M3"),
  data.frame(N = model_outputs[[4]][["PosteriorAllRegions"]], model = "M4"),
  data.frame(N = model_outputs[[5]][["PosteriorAllRegions"]], model = "M5"),
  data.frame(N = model_outputs[[6]][["PosteriorAllRegions"]], model = "SCR")
)

N_original <- 952

res_bias_ntot_plot <- N_all %>%
  mutate(model = factor(model, levels = c("M1", "M2", "M3", "M4", "M5", "SCR"))) %>%
  ggplot(aes(x = model, y = N, fill = model)) +
  # SCR reference lines
  geom_hline(yintercept = N_original,alpha = 0.9, linetype ="dashed", color = "#D7263D", size = 1) +
  geom_violin(color = NA, trim = FALSE) +
  stat_summary(fun = median, geom = "point", color = "white", size = 0.5) +
  scale_y_continuous(limits = c(500, 3500), breaks = seq(500, 3500, 500)) +
  scale_fill_manual(values = c(
    "M1" = "#ccd5dd",
    "M2" = "#EBCC2A",
    "M3" = "#E1AF00",
    "M4" = "skyblue",
    "M5" = "#3B9AB2",
    "SCR" = "#D7263D"
  )) +
  labs(y = "Population size", x = "") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.5, 'cm'),
    legend.key.width = unit(0.5, 'cm'),
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 10)
  )

# Print it
res_bias_ntot_plot
ggsave(
  filename = "/Users/virginia/Desktop/Wolf_Npost_Violin.png",
  plot = res_bias_ntot_plot,
  width = 8,           # inches
  height = 5,        # adjust as needed
  dpi = 300,
  units = "in"
)




