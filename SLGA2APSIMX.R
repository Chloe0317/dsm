# File: DSM2APSIMX.R
# Original Author: Chloe Lai (Chloe.Lai@unisq.edu.au)
# License: Creative Commons Attribution 4.0 International (CC BY 4.0)
# URL: https://creativecommons.org/licenses/by/4.0/
# Description: This script performs extraction of soil information (e.g., sand, silt, clay content)
# based on coordinates or a shapefile from the Soil and Landscape Grid of Australia and write it into apsimx


# to install ithir, uncomment the following
# install.packages("devtools")
# 
# library(devtools)
# 
# install_bitbucket("brendo1001/ithir/pkg")

suppressPackageStartupMessages({
  library(ithir)
  library(terra)
  library(apsimx)
})

# functions ---------
get_slga_property <- function(property, namestring, version, depthlayer = 1:6, extent = NULL, apikey) { 
  # Define depth layers (predefined soil depths in cm)
  depths <- c("_000_005_", "_005_015_", "_015_030_", "_030_060_", "_060_100_", "_100_200_")
  # Select specified depth layers
  depths <- depths[depthlayer]
  
  # Construct base URL for SLGA data with API key authentication
  url1 <- paste0("/vsicurl/https://", apikey, "@data.tern.org.au/model-derived/slga/NationalMaps/SoilAndLandscapeGrid/")
  
  # Append property, version, and depth layers to create full URL paths
  if (property == "pHc") {
    url2 <- paste0(url1, property, "/v", version, "/", "PHC", depths, namestring)
  } else if (property == "SOC") {
    url2 <- paste0(url1, property, "/90m/v", version, "/", "SOC", depths, namestring)
  } else {
    url2 <- paste0(url1, property, "/v", version, "/", property, depths, namestring)
  }
  
  # Retrieve and crop each depth layer if extent is provided
  rast_list <- lapply(url2, function(url) {
    r <- rast(url)  # Load raster
    if (!is.null(extent)) {
      # Crop to the specified extent
      r <- crop(r, extent)
    }
    return(r)
  })
  
  # Stack the individual rasters into a single raster stack
  rastStack <- do.call(c, rast_list)
  
  return(rastStack)
}




# SLGA -----

####  Access a COG requiring authentication from the TERN Data Store  
# please create your own API from Tern, information can be found here
# https://esoil.io/TERNLandscapes/Public/Pages/SLGA/GetData-COGSDataStore.html

apistring <- "xxxyourternapikeyxxxx"
apikey <- paste0('apikey:', apistring)


# extract at point location soil information using the coordinates of Dalby in the Wheat.apsimx Example 
xy <- c(151.26, -27.18)
slga <- read.table("SLGA.csv", sep = ",", header = T)
# Initialize property_df with 6 rows (depths) and 7 columns (properties)
property_df <- data.frame(matrix(NA, nrow = 6, ncol = 7))
colnames(property_df) <- slga$SLGAProperty  # Set column names as property names

for (i in (1:nrow(slga))){
  rr <- get_slga_property(slga$SLGAProperty[i], slga$URL[i], slga$Version[i], apikey = apikey)
  property.i <- extract(rr, xy)
  property_df[[slga$SLGAProperty[i]]] <- property.i
  
}

# for field-based extraction, calculate the average value across the polygon
myfield <- vect("myfield.shp")
# Initialize property_df with 6 rows (depths) and 7 columns (properties)
property_df <- data.frame(matrix(NA, nrow = 6, ncol = 7))
colnames(property_df) <- slga$SLGAProperty  # Set column names as property names

for (i in 1:nrow(slga)) {
  # Retrieve SLGA property raster
  rr <- get_slga_property(slga$SLGAProperty[i], slga$URL[i], slga$Version[i], extent = ext(myfield), apikey = apikey)
  
  # Calculate layer means using terra::extract
  layer_means <- terra::extract(rr, vect(burdekin), fun = "mean", exact = TRUE)
  
  # Process and round the extracted means
  if (slga$SLGAProperty[i] %in% c("DUL", "L15")) {
    values <- round(as.numeric(layer_means)[2:7] / 100, digits = 3)
  } else {
    values <- round(as.numeric(layer_means)[2:7], digits = 2)
  }
  
  # Add the results as a new column in the data frame
  property_df[[slga$SLGAProperty[i]]] <- values
}

# harmonise to the required depth interval for your APSIM soil profile 
# for demonstration here I used the Wheat.apsimx in the example folder


# ea_spline requires the first three columns to be ID, UpperDepth and LowerDepth
# Define GlobalSoilMap soil depth intervals
depth_intervals <- c("0-5 cm", "5-15 cm", "15-30 cm", "30-60 cm", "60-100 cm", "100-200 cm")

# Split the intervals into upper and lower depths
depths <- do.call(rbind, strsplit(gsub(" cm", "", depth_intervals), "-"))

# Create the data frame
depth_df <- data.frame(
  ID = 1,
  UpperDepth = as.numeric(depths[, 1]), # Upper depth
  LowerDepth = as.numeric(depths[, 2])  # Lower depth
)

property_df <- cbind.data.frame(depth_df, property_df)

# the Wheat example uses APsoil 104 that have the following depth layers, 
# 0-15, 15-30, 30-60, 60-90, 90-120, 120-150, 150-180
# Initialize harmonised_df with 7 rows (depths) and 7 columns (properties)
depths <- c(0, 15, 30, 60, 90, 120, 150, 180) 
harmonised_df <- data.frame(matrix(NA, nrow = 7, ncol = 7))
colnames(harmonised_df) <- slga$SLGAProperty  # Set column names as property names

for (i in (slga$SLGAProperty)){

  harmonised_values <- as.numeric(ea_spline(property_df, i, d = depths, show.progress = F)$harmonised[,2:(length(depths))])
  # Add the results as a new column in the data frame
  harmonised_df[[i]] <- round(harmonised_values, digits = 2)
}

# or directly write into your apsim soil profile 
# (make sure you have a soil set up with the GlobalSoilMap depth intervals)
# as SLGA currently doesn't contain SAT you might want to check to ensure your SAT > DUL
# you can also use the apsimx_soil_profile function to create a soil profile in R

#Specify the apsimx parameter and path

param_paths <- c(
  "BD",
  "ParticleSizeClay",
  "ParticleSizeSilt",
  "ParticleSizeSand",
  "DUL",
  "LL15", 
  "PH"
)

soil_physical <- c(
  "BD",
  "ParticleSizeClay",
  "ParticleSizeSilt",
  "ParticleSizeSand",
  "DUL",
  "LL15"
)
soil_chemical <- c("PH")

# Loop through each soil property and update APSIM
for (i in (1:length(param_paths))) {
  
  if (param_paths[i]%in%soil_physical){
    edit_apsimx(
      file = "Wheat.apsimx",     
      node = "Soil",                 
      soil.child = "Physical",         
      parm = param_paths[i],  
      value = harmonised_df[, i],      
      verbose = TRUE,
      overwrite = TRUE)
  }else{
    edit_apsimx(
      file = "Wheat.apsimx",         
      node = "Soil",                 
      soil.child = "Chemical",         
      parm = param_paths[i],  
      value = harmonised_df[, i],      
      verbose = TRUE,
      overwrite = TRUE)
  }
  
}
