# File: SoilTextureExtraction.R
# Original Author: Chloe Lai (Chloe.Lai@unisq.edu.au)
# License: Creative Commons Attribution 4.0 International (CC BY 4.0)
# URL: https://creativecommons.org/licenses/by/4.0/
# Description: This script performs extraction of texture information (sand, silt, clay content)
# based on coordinates from the Soil and Landscape Grid of Australia and from the SoilGrids.

library(terra)
# functions ---------
get_slga_texture <- function(property, apikey){
  # function to get the slga raster sets 
  # the namestring is the last bit of naming convention for the texture maps
  namestring <- "EV_N_P_AU_TRN_N_20210902.tif"
  depths <- c("_000_005_", "_005_015_", "_015_030_", "_030_060_", "_060_100_", "_100_200_")
  url1 <- paste0("/vsicurl/https://", apikey, "@data.tern.org.au/model-derived/slga/NationalMaps/SoilAndLandscapeGrid/")
  url2 <- paste0(url1, property, "/v2/", property, depths, namestring)
  rast_list <- lapply(url2, rast)
  rastStack <- do.call(c, rast_list)
}

get_soilgrids_texture <- function(voi, 
                                  soilDepth = c('0-5cm', "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm"),
                                  minlong, maxlong, minlat, maxlat){
  url11 <- "https://maps.isric.org/mapserv?map=/map/"
  url12 <- ".map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID="
  extractionString <- paste(voi, soilDepth, "mean", sep = "_")
  url21 <- "&FORMAT=image/tiff&SUBSET=long("
  bbox <- paste0(minlong, ",", maxlong, ")&SUBSET=lat(", minlat, ",", maxlat)
  url22 <- paste0(bbox, ")&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326")
  url <- paste0(url11, voi, url12, extractionString, url21, url22)
  rast_list <- lapply(url, rast)
  rastStack <- do.call(c, rast_list)
  
}


# SLGA -----

####  Access a COG requiring authentication from the TERN Data Store  
# please create your own API from Tern, information can be found here
# https://esoil.io/TERNLandscapes/Public/Pages/SLGA/GetData-COGSDataStore.html

apistring <- "xxxyourternapikeyxxxx"
apikey <- paste0('apikey:', apistring)


# soil depths of SLGA based on GlobalSoilMap Project, 
# don't use space as we'll use it to call SoilGrids data as well
soilDepth <- c('0-5cm', "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm")

# for texture, SLGA designation is CLY, SND and SLT for clay, sand and silt and 
# the naming convention to use is otherwise the same

# Create a dataframe with coordinates for Australian state capitals
state_capitals <- data.frame(
  City = c("Sydney", "Melbourne", "Brisbane", "Perth", "Adelaide", "Hobart", "Darwin", "Canberra"),
  State = c("New South Wales", "Victoria", "Queensland", "Western Australia", "South Australia", "Tasmania", "Northern Territory", "Australian Capital Territory"),
  Longitude = c(151.2093, 144.9631, 153.0251, 115.8605, 138.6007, 147.3272, 130.8456, 149.1300),
  Latitude = c(-33.8688, -37.8136, -27.4698, -31.9505, -34.9285, -42.8821, -12.4634, -35.2809)
  
)

clay.slga.r <- get_slga_texture("CLY", apikey)
sand.slga.r <- get_slga_texture("SND", apikey)
silt.slga.r <- get_slga_texture("SLT", apikey)
clay_slga <- extract(clay.slga.r, state_capitals[,3:4])
sand_slga <- extract(sand.slga.r, state_capitals[,3:4])
silt_slga <- extract(silt.slga.r, state_capitals[,3:4])
colnames(clay_slga)[2:7] <- colnames(sand_slga)[2:7] <- colnames(silt_slga)[2:7] <- soilDepth

# WoSIS  SoilGrids----
# url https://soilgrids.org/
# This can also be accessed via GEE 
# further info https://git.wur.nl/isric/soilgrids/soilgrids.notebooks/-/blob/master/markdown/access_on_gee.md
# variable naming convention can be found https://www.isric.org/explore/soilgrids/faq-soilgrids#What_do_the_filename_codes_mean
# for texture, they are just "clay, sand, silt"


# constrain it to a bounding box for illustration 
# you can change the bounding box based on points you'd like to extract
clay <- get_soilgrids_texture(voi = "clay", minlong = 152, maxlong = 153, minlat = -27, maxlat = -26)
sand <- get_soilgrids_texture(voi = "sand", minlong = 152, maxlong = 153, minlat = -27, maxlat = -26)
silt <- get_soilgrids_texture(voi = "silt", minlong = 152, maxlong = 153, minlat = -27, maxlat = -26)

# note that soilgrids textures are reported as g/kg, so convert to %, divide by 10

clay_soilgrids <- extract(clay, data.frame(152.05, -26.05))/10
sand_soilgrids <- extract(sand, data.frame(152.05, -26.05))/10
silt_soilgrids <- extract(silt, data.frame(152.05, -26.05))/10
colnames(clay_soilgrids)[2:7] <- colnames(sand_soilgrids)[2:7] <- colnames(silt_soilgrids)[2:7] <- soilDepth

