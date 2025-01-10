library(httr)
library(jsonlite)
library(apsimx)
# Define the function
get_slga_soil <- function(latitude, longitude) {
  # Define properties
  properties <- list(
    clay = "https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_000_005_EV_N_P_AU_TRN_N_20210902.tif",
    sand = "https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/SND/SND_000_005_EV_N_P_AU_TRN_N_20210902.tif",
    silt = "https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/SLT/SLT_000_005_EV_N_P_AU_TRN_N_20210902.tif",
    wv1500 = "https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/L15/L15_000_005_EV_N_P_AU_TRN_N_20210614.tif",
    wv0033 = "https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/DUL/DUL_000_005_EV_N_P_AU_TRN_N_20210614.tif",
    bdod = "https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/BDW/BDW_000_005_EV_N_P_AU_TRN_N_20230607.tif",
    nitrogen = "https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/NTO/NTO_000_005_EV_N_P_AU_NAT_C_20231101.tif",
    phh2o = "https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/PHW/PHW_000_005_EV_N_P_AU_TRN_N_20220520.tif",
    cec = "https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CEC/CEC_000_005_EV_N_P_AU_TRN_N_20220826.tif",
    des = "https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/DES/DES_000_200_EV_N_P_AU_TRN_C_20190901.tif",
    soc = "https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/SOC/SOC_000_005_EV_N_P_AU_TRN_N_20220727.tif"
  )
  
  # Define depths
  depths <- c("000_005", "005_015", "015_030", "030_060", "060_100", "100_200")
  
  # Initialize results list
  results <- list()
  
  # Loop through each property and depth
  for (property_name in names(properties)) {
    property_url <- properties[[property_name]]
    property_results <- numeric(length(depths))  # Initialize vector of length 6
    
    for (i in seq_along(depths)) {
      depth <- depths[i]
      
      # Replace "000_005" with the current depth
      cog_path <- gsub("000_005", depth, property_url)
      
      # Construct the API request URL
      api_url <- paste0(
        "https://esoil.io/TERNLandscapes/RasterProductsAPI/Drill?",
        "format=json&verbose=false&COGPath=", cog_path,
        "&latitude=", latitude, "&longitude=", longitude
      )
      
      # Make the API request
      response <- GET(api_url)
      
      # Parse the response
      if (status_code(response) == 200) {
        data <- fromJSON(content(response, "text"))
        property_results[i] <- data$Value  # Assign to vector
      } else {
        property_results[i] <- NA  # Handle errors gracefully
      }
    }
    
    # Store the results for the property
    results[[property_name]] <- property_results
  }
  
  # Return the results
  return(results)
}



slga <- get_slga_soil(lat, lon)
thcknss <- c(50, 100, 150, 300, 400, 1000) ## in mm

soil_profile <- apsimx_soil_profile(nlayers = 6, Thickness = thcknss) 
soil_profile$soil$ParticleSizeClay <- NA
soil_profile$soil$ParticleSizeSilt <- NA
soil_profile$soil$ParticleSizeSand <- NA
soil_profile$soil$CEC <- NA  
soil_profile$soil$Nitrogen <- NA
soil_profile$soil$DUL <- NA
soil_profile$soil$LL15 <- NA
soil_profile$soil$SAT <- NA

# Assigning the entire vector from slga to soil_profile
soil_profile$soil$BD <- slga[["bdod"]]  # Bulk density
soil_profile$soil$Carbon <- slga[["soc"]]  # Soil organic carbon
soil_profile$soil$PH <- slga[["phh2o"]]  # pH
soil_profile$soil$ParticleSizeClay <- slga[["clay"]]  # Clay percentage
soil_profile$soil$ParticleSizeSand <- slga[["sand"]]  # Sand percentage
soil_profile$soil$ParticleSizeSilt <- slga[["silt"]] # Silt percentage
soil_profile$soil$Nitrogen <- slga[["nitrogen"]]  # Nitrogen content
soil_profile$soil$CEC <- slga[["cec"]]  # Cation exchange capacity
soil_profile$soil$DUL <- slga[["wv0033"]]* 1e-2  # Drained upper limit
soil_profile$soil$LL15 <- slga[["wv1500"]]* 1e-2  # Lower limit (wilting point)

# Calculate SAT (saturation) for each depth layer based on BD and particle density (2.65 g/cm³)
particle_density <- 2.65  # Particle density (g/cm³)
soil_profile$soil$SAT <- 1 - (soil_profile$soil$BD / particle_density)

# Applying pedotransfer function for KS in Saxton and Rawls, 2006.
B <- (log(1500) - log(33))/(log(soil_profile$soil$DUL) - log(soil_profile$soil$LL15))
Lambda <- 1/B
soil_profile$soil$KS <- (1930 * (soil_profile$soil$SAT - soil_profile$soil$DUL)^(3 - Lambda)) * 100

soil_profile$soil$AirDry <- soil_profile$soil$LL15
soil_profile$soil$AirDry[1] <- soil_profile$soil$LL15[1] * 0.5 ## AirDry is half of LL for the first layer

for(i in soil_profile$crops){
  soil_profile$soil[[paste0(i,".LL")]] <- soil_profile$soil$LL15 ## Without better information  
}

#### Passing parameters from soilwat
## The soil texture class in the metadata will be based on the first layer only
# note that this is the USA texture triangle - will need to update to the Australian one later
txt_clss <- texture_class(soil_profile$soil$ParticleSizeClay * 1e-2, soil_profile$soil$ParticleSizeSilt * 1e-2)


# redefine the texture_class function to accept vectors
texture2soilParms <- function(texture.class = "NO DATA") { 
  # Define texture classes and associated parameters
  textureClasses <- c("clay", "silty clay", "sandy clay", "clay loam", "silty clay loam", 
                      "sandy clay loam", "loam", "silty loam", "sandy loam", "silt", 
                      "loamy sand", "sand", "NO DATA")  
  Albedo <- c(0.12, 0.12, 0.13, 0.13, 0.12, 0.13, 0.13, 0.14, 0.13, 0.13, 0.16, 0.19, 0.13)
  CN2 <- c(73.0, 73.0, 73.0, 73.0, 73.0, 73.0, 73.0, 73.0, 68.0, 73.0, 68.0, 68.0, 73.0)
  SWCON <- c(0.25, 0.3, 0.3, 0.4, 0.5, 0.5, 0.5, 0.5, 0.6, 0.5, 0.6, 0.75, 0.5)
  
  # Match the first layer of texture.class
  wtc <- match(texture.class[1], textureClasses, nomatch = length(textureClasses))
  
  # Match all layers of texture.class for SWCON
  swcon <- SWCON[match(texture.class, textureClasses, nomatch = length(textureClasses))]
  
  # Return results
  ans <- list(
    textureClasses = textureClasses[wtc],  # Only for the first layer
    Albedo = Albedo[wtc],                  # Only for the first layer
    CN2 = CN2[wtc],                        # Only for the first layer
    SWCON = swcon                          # Vector for all layers
  )
  return(ans)
}


# Get soil parameters
t2sp <- texture2soilParms(txt_clss)

# get the Australian Soil Order 
# it can be accessed here, however, I don't think point-drill is allowed yet
# https://swift.rc.nectar.org.au/v1/AUTH_05bca33fce34447ba7033b9305947f11/landscapes-csiro-slga-public/NationalMaps/SoilClassifications/ASC/90m/ASC_EV_C_P_AU_TRN_N.cog.tif

#### Attributes ####
alist <- list()
alist$SoilType <- paste("SoilType = ", txt_clss)
alist$State <- state
alist$Country <- country
alist$Longitude <- lon
alist$Latitude <- lat
alist$DataSource <- paste("Original source is Soil and Landscape Grid of Australia www.isric.org. See: https://esoil.io/TERNLandscapes/Public/Pages/SLGA/index.html",Sys.time())
alist$Comments <- paste("resolution = 90 m",
                        "- taxonomic classification name =", txt_clss,
                        "- drainage class =", NA, 
                        "- elevation =", NA,
                        "- slope =", NA,
                        "- geomdesc =", NA)

soil_profile$metadata <- alist