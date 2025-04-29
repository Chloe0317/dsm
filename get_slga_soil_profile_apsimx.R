# this script is to download an apsimx-compatible SLGA soil profile 
# author: Chloe Lai
# some functions are taken/adapted from the r apsimx package written by Fernando Miguez 
# some are from Andrew Moore (World Modeller Database)
# there are inconsistencies among the approaches to approx. different soil properties 
# and min. efforts have yet to be made to consolidate these

# the soil specific crop parameters (CLL, KL and XF) are calculated using Andrew's method
# with an additional condition where a soil with high clay content the rooting depth is reduced by 300 mm
# as per Dagliesh et al. (2016)
# (defined as clay% > 40 throughout of the soil profile)

suppressPackageStartupMessages({
  library(apsimx)
  library(jsonlite)
  library(curl)
  library(httr)
  
})

# required functions -------------------------
# functions from Andrew Moore
# Re-express the PSD in terms of the International system, using an equation from Minasny et al. (2001)

intl_clay_propn <- function( usda_clay, usda_silt ) { 
  return( usda_clay) 
}

intl_silt_propn <- function( usda_clay, usda_silt ) { 
  return( max( 0.0, -0.0041 - 0.127*usda_clay + 0.553*usda_silt + 0.17*usda_clay^2 - 0.19*usda_silt^2 + 0.59*usda_clay*usda_silt ) ) 
}

intl_sand_propn <- function( usda_clay, usda_silt ) { 
  return( 1.0 - intl_clay_propn( usda_clay, usda_silt ) - intl_silt_propn( usda_clay, usda_silt ) )
}  


## Texture to other parameters

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

#==========================================================================
# Calculation of soil-crop parameters
#==========================================================================

# Parameters for the XF equation


crop_types         <- c("Rice","Wheat","Sugarcane","Maize","Soybean","OilPalm","Cassava","Plantain","Sunflower","Canola","Sorghum","Millet","Peanut","Barley","Cotton", "Oats", "Chickpea", "Mungbean", "Pigeon Pea", "Cowpea", "Fababean", "Lucerne", "Lupin", "Lentil", "Triticale", "SCRUM")
kl_max             <- c(  0.06,   0.06,       0.06,   0.06,     0.06,     0.06,     0.06,      0.06,       0.10,    0.06,     0.07,    0.07,    0.06,    0.07,    0.10,   0.06,       0.06,       0.06,         0.06,     0.06,       0.08,      0.10,    0.06,     0.06,        0.07, 0.06 )
depth_max_kl       <- c(   400,    400,        400,    400,      400,      400,      400,       400,        100,     400,      400,     400,     600,     400,     800,    400,        600,        500,          600,      400,        700,      1300,     400,      400,         300, 400)
depth_max_root     <- c(  1900,   1900,       3000,   1900,     1900,     3000,     1900,      3000,       1800,    1900,     2400,    2400,    1800,    2100,    2600,   1900,       1500,       1500,         1800,     1900,       2100,      3500,    1900,     1500,        1400, 1900)

bd_restrict_sand_0 <- rep( 1.40, length(crop_types) )  # BD below which XF=1 when propn sand = 0.0
bd_restrict_sand_1 <- rep( 1.60, length(crop_types) )  # BD below which XF=1 when propn sand = 1.0
bd_zero_sand_0     <- rep( 1.80, length(crop_types) )  # BD above which XF=0 when propn sand = 0.0 (notionally)
crop_min_xf        <- rep( 0.10, length(crop_types) )

xf_param_1         <- bd_restrict_sand_0
xf_param_2         <- bd_restrict_sand_1 - bd_restrict_sand_0
xf_param_3         <- bd_zero_sand_0     - bd_restrict_sand_0

compute_crop_params <- function( crop, thickness, bd_fine, intl_sand_propn, intl_clay_propn, ll15, dul ) {
  
  if (all(intl_clay_propn > 0.4)){
    depth_max_root <- depth_max_root - 300
  }
  
  layer_count     <- length(thickness)
  crop_id         <- match( crop, crop_types)
  crop_kl_max     <- kl_max[crop_id]  
  crop_depth_kl   <- depth_max_kl[crop_id]
  crop_depth_root <- depth_max_root[crop_id]
  p1              <- xf_param_1[crop_id]
  p2              <- xf_param_2[crop_id]
  p3              <- xf_param_3[crop_id]
  p4              <- crop_min_xf[crop_id]
  
  # Andrew Moore's equation for xf
  xf <- 1.0 - (bd_fine - (p1 + p2 * intl_sand_propn)) / p3
  xf <- ifelse( xf<p4, p4, ifelse( xf<1.0, xf, 1.0 ) )
  
  kl             <- rep( 0.0, layer_count )
  ll             <- rep( 0.0, layer_count )
  top_eff_depth  <- 0.0
  for (layer in 1:layer_count) {
    eff_thickness      <- thickness[layer] / xf[layer]
    bottom_eff_depth   <- top_eff_depth + eff_thickness
    
    propn_at_max_kl    <- max(      0.0, min( bottom_eff_depth, crop_depth_kl   ) - top_eff_depth ) / eff_thickness
    propn_at_decr_kl   <- max( max( 0.0, min( bottom_eff_depth, crop_depth_root ) - top_eff_depth ) / eff_thickness - propn_at_max_kl, 0.0 )
    propn_at_zero_kl   <- 1.0 - propn_at_max_kl - propn_at_decr_kl
    
    ratio_at_top_depth    <- max( 0.0, min( (crop_depth_root-top_eff_depth   )/(crop_depth_root-crop_depth_kl), 1.0 ) )
    ratio_at_bottom_depth <- max( 0.0, min( (crop_depth_root-bottom_eff_depth)/(crop_depth_root-crop_depth_kl), 1.0 ) )
    mean_decreasing_ratio <- 0.5 * (ratio_at_top_depth + ratio_at_bottom_depth)
    
    weighted_ratio   <- propn_at_max_kl * 1.0 + propn_at_decr_kl * mean_decreasing_ratio + propn_at_zero_kl * 0.0
    
    kl[layer]        <- crop_kl_max * weighted_ratio
    ll[layer]        <- ll15[layer] + (dul[layer] - ll15[layer]) * (1.0 - weighted_ratio)
    
    top_eff_depth    <- bottom_eff_depth
  }
  
  # Set xf=0.0 below the rooting zone for clarity
  xf <- xf * (kl > 0.0)  
  
  # Create and return a dataframe
  result <- data.frame(
    setNames(
      list(kl, ll, xf),
      paste0(crop, c(".KL", ".LL", ".XF"))
    )
  )
  return(result)
}

compute_all_crops_params <- function(crops, thickness, bd_fine, intl_sand_propn, intl_clay_propn, ll15, dul) {
  result_list <- lapply(crops, function(crop) {
    compute_crop_params(crop, thickness, bd_fine, intl_sand_propn, intl_clay_propn, ll15, dul)
  })
  
  result_df <- do.call(cbind, result_list)
  return(result_df)
}

approx_soil_variable <- function(x, xout = NULL, soil.bottom = 200, 
                                 method = c("constant", "linear"),
                                 nlayers = 10){
  
  xd <- as.data.frame(x)
  method <- match.arg(method)
  if(is.null(xout)) xout <- seq(0, soil.bottom, length.out = nlayers)
  
  ans <- stats::approx(x = xd[[1]], y = xd[[2]], method = method, xout = xout, rule = 2)
  
  ansd <- data.frame(x = ans$x, y = ans$y)
  return(ansd)
}

# redefine the texture_class function to accept vectors
texture_class <- function(usda_clay, usda_silt) {
  if (any(usda_clay < 0 | usda_clay > 1)) stop("All values in usda_clay should be between 0 and 1")
  if (any(usda_silt < 0 | usda_silt > 1)) stop("All values in usda_silt should be between 0 and 1")
  
  intl_clay <- usda_clay
  intl_silt <- usda_silt
  intl_sand <- 1.0 - intl_clay - intl_silt
  
  # Initialize result vector
  classes <- rep(NA, length(usda_clay))
  
  # Apply texture triangle rules
  classes[(intl_sand < 0.75 - intl_clay) & (intl_clay >= 0.40)] <- "silty clay"
  classes[(intl_sand < 0.75 - intl_clay) & (intl_clay >= 0.26) & is.na(classes)] <- "silty clay loam"
  classes[(intl_sand < 0.75 - intl_clay) & is.na(classes)] <- "silty loam"
  classes[(intl_clay >= 0.40 + (0.305 - 0.40) / (0.635 - 0.35) * (intl_sand - 0.35)) & 
            (intl_clay < 0.50 + (0.305 - 0.50) / (0.635 - 0.50) * (intl_sand - 0.50)) & is.na(classes)] <- "clay"
  classes[(intl_clay >= 0.26 + (0.305 - 0.26) / (0.635 - 0.74) * (intl_sand - 0.74)) & is.na(classes)] <- "sandy clay"
  classes[(intl_clay >= 0.26 + (0.17 - 0.26) / (0.83 - 0.49) * (intl_sand - 0.49)) & 
            (intl_clay < 0.10 + (0.305 - 0.10) / (0.635 - 0.775) * (intl_sand - 0.775)) & is.na(classes)] <- "clay loam"
  classes[(intl_clay >= 0.26 + (0.17 - 0.26) / (0.83 - 0.49) * (intl_sand - 0.49)) & is.na(classes)] <- "sandy clay loam"
  classes[(intl_clay >= 0.10 + (0.12 - 0.10) / (0.63 - 0.775) * (intl_sand - 0.775)) & 
            (intl_clay < 0.10 + (0.305 - 0.10) / (0.635 - 0.775) * (intl_sand - 0.775)) & is.na(classes)] <- "loam"
  classes[(intl_clay >= 0.10 + (0.12 - 0.10) / (0.63 - 0.775) * (intl_sand - 0.775)) & is.na(classes)] <- "sandy loam"
  classes[(intl_clay < 0.00 + (0.08 - 0.00) / (0.88 - 0.93) * (intl_sand - 0.93)) & is.na(classes)] <- "loamy sand"
  classes[is.na(classes)] <- "sand"
  
  return(classes)
}

calculate_SAT <- function(texture.class, BD) {
  
  # Define the adjustment factor (e) based on broad soil categories
  e_values <- c("clay" = 0.03, 
                "silty clay" = 0.03, 
                "sandy clay" = 0.03, 
                "clay loam" = 0.05, 
                "silty clay loam" = 0.05, 
                "sandy clay loam" = 0.05, 
                "loam" = 0.05, 
                "silty loam" = 0.05, 
                "sandy loam" = 0.05, 
                "loamy sand" = 0.07, 
                "sand" = 0.07)
  
  # Map the adjustment factor (e) based on the texture.class
  e <- e_values[texture.class]
  
  # Calculate Porosity (PO) as a percentage
  PO <- 1 - BD / 2.65
  
  # Calculate Saturation (SAT) as a percentage
  SAT <- PO - e 
  
  # Return SAT values
  return(SAT)
}

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
    # des = "https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/DES/DES_000_200_EV_N_P_AU_TRN_C_20190901.tif",
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
        data <- suppressMessages(fromJSON(content(response, "text")))
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

get_slga_soil_profile <- function(lonlat, 
                                  statistic = c("mean", "Q0.5"),
                                  soil.profile,
                                  find.location.name = TRUE,
                                  fix = FALSE,
                                  verbose = TRUE,
                                  check = TRUE,
                                  physical = c("default", "SR"),
                                  xargs = NULL){
  
  statistic <- match.arg(statistic)
  physical <- match.arg(physical)
  
  #### Create extent step ####
  lon <- as.numeric(lonlat[1])
  lat <- as.numeric(lonlat[2])
  
  if (lon < 112 || lon > 154) stop("Longitude should be between 112 and 154 for the extent of Australia")
  if (lat < -44 || lat > -10) stop("Latitude should be between -44 and -10 for the extent of Australia")
  
  slga <- get_slga_soil(lat, lon)
  
  ## These are the default thicknesses in SLGA
  thcknss <- c(50, 100, 150, 300, 400, 1000) ## in mm
  
  ## Some variables can be passed to apsimx:::approx_soil_variable
  soil.bottom <- 200
  method <- "constant"
  nlayers <- 10
  crps <- c("Maize", "Soybean", "Wheat")
  if(!is.null(xargs)){
    ### Soil bottom
    if(!is.null(xargs$soil.bottom)){
      soil.bottom <- xargs$soil.bottom
    }
    ### Method
    if(!is.null(xargs$method)){
      method <- xargs$method
    }
    ### Number of layers
    if(!is.null(xargs$nlayers)){
      nlayers <- xargs$nlayers
    }
    ### Crops
    if(!is.null(xargs$crops)){
      crps <- xargs$crops
    }
  }
  
  ## Create the empty soil profile
  if(missing(soil.profile)){
    new.soil <- FALSE
    soil_profile <- apsimx_soil_profile(nlayers = 6, Thickness = thcknss, crops = crps) 
    soil_profile$soil$ParticleSizeClay <- NA
    soil_profile$soil$ParticleSizeSilt <- NA
    soil_profile$soil$ParticleSizeSand <- NA
    soil_profile$soil$CEC <- NA  
    soil_profile$soil$Nitrogen <- NA
    soil_profile$soil$DUL <- NA
    soil_profile$soil$LL15 <- NA
    soil_profile$soil$SAT <- NA
    
  }else{
    ## stop("This is not fully implemented yet. Submit a github issue if you need it.", call. = FALSE)
    soil_profile <- soil.profile
    new.soil <- TRUE
  }
  
  
  
  if (new.soil) {
    sp.xout <- cumsum(soil_profile$soil$Thickness)
    
    soil_profile$soil$BD <- approx_soil_variable(data.frame(x = cumsum(thcknss), y = slga[["bdod"]]), 
                                                 xout = sp.xout, soil.bottom = soil.bottom, method = method, nlayers = nlayers)$y    
    soil_profile$soil$Carbon <- approx_soil_variable(data.frame(x = cumsum(thcknss), y = slga[["soc"]]), 
                                                     xout = sp.xout, soil.bottom = soil.bottom, method = method, nlayers = nlayers)$y    
    soil_profile$soil$PH <- approx_soil_variable(data.frame(x = cumsum(thcknss), y = slga[["phh2o"]]), 
                                                 xout = sp.xout, soil.bottom = soil.bottom, method = method, nlayers = nlayers)$y    
    soil_profile$soil$ParticleSizeClay <- approx_soil_variable(data.frame(x = cumsum(thcknss), y = slga[["clay"]]), 
                                                               xout = sp.xout, soil.bottom = soil.bottom, method = method, nlayers = nlayers)$y   
    soil_profile$soil$ParticleSizeSand <- approx_soil_variable(data.frame(x = cumsum(thcknss), y = slga[["sand"]]), 
                                                               xout = sp.xout, soil.bottom = soil.bottom, method = method, nlayers = nlayers)$y 
    soil_profile$soil$ParticleSizeSilt <- approx_soil_variable(data.frame(x = cumsum(thcknss), y = slga[["silt"]]), 
                                                               xout = sp.xout, soil.bottom = soil.bottom, method = method, nlayers = nlayers)$y 
    soil_profile$soil$Nitrogen <- approx_soil_variable(data.frame(x = cumsum(thcknss), y = slga[["nitrogen"]]), 
                                                       xout = sp.xout, soil.bottom = soil.bottom, method = method, nlayers = nlayers)$y   
    soil_profile$soil$CEC <- approx_soil_variable(data.frame(x = cumsum(thcknss), y = slga[["cec"]]), 
                                                  xout = sp.xout, soil.bottom = soil.bottom, method = method, nlayers = nlayers)$y 
    
    soil_profile$soil$DUL <- approx_soil_variable(data.frame(x = cumsum(thcknss), y = slga[["wv0033"]] ), 
                                                  xout = sp.xout, soil.bottom = soil.bottom, method = method, nlayers = nlayers)$y * 1e-2  
    soil_profile$soil$LL15 <- approx_soil_variable(data.frame(x = cumsum(thcknss), y = slga[["wv1500"]]), 
                                                   xout = sp.xout, soil.bottom = soil.bottom, method = method, nlayers = nlayers)$y * 1e-2   
  } else {
    soil_profile$soil$BD <- slga[["bdod"]]  # Bulk density
    soil_profile$soil$Carbon <- slga[["soc"]]  # Soil organic carbon
    soil_profile$soil$PH <- slga[["phh2o"]]  # pH
    soil_profile$soil$ParticleSizeClay <- slga[["clay"]]  # Clay percentage
    soil_profile$soil$ParticleSizeSand <- slga[["sand"]]  # Sand percentage
    soil_profile$soil$ParticleSizeSilt <- slga[["silt"]]  # Silt percentage
    soil_profile$soil$Nitrogen <- slga[["nitrogen"]]  # Nitrogen content
    soil_profile$soil$CEC <- slga[["cec"]]  # Cation exchange capacity
    soil_profile$soil$DUL <- slga[["wv0033"]] * 1e-2  # Drained upper limit
    soil_profile$soil$LL15 <- slga[["wv1500"]] * 1e-2  # Lower limit (wilting point)
    
    
  }
  
  
  
  ## Populating DUL and LL. These are equations from Table 1 in Saxton and Rawls 2006
  if(physical == "SR"){
    soil_profile$soil$DUL <- sr_dul(soil_profile$soil$ParticleSizeClay, soil_profile$soil$ParticleSizeSand, soil_profile$soil$Carbon * 2)
    soil_profile$soil$LL15 <- sr_ll(soil_profile$soil$ParticleSizeClay, soil_profile$soil$ParticleSizeSand, soil_profile$soil$Carbon * 2)
    DUL_S <- sr_dul_s(soil_profile$soil$ParticleSizeClay, soil_profile$soil$ParticleSizeSand, soil_profile$soil$Carbon * 2)
    soil_profile$soil$SAT <- sr_sat(soil_profile$soil$ParticleSizeSand, soil_profile$soil$DUL, DUL_S)    
  }
  
  
  soil_profile$soil$AirDry <- soil_profile$soil$LL15
  soil_profile$soil$AirDry[1] <- soil_profile$soil$LL15[1] * 0.5 ## AirDry is half of LL for the first layer
  soil_profile$soil$AirDry[2] <- soil_profile$soil$LL15[2] * 0.8 ## AirDry is 80% of LL until 60 cm
  soil_profile$soil$AirDry[3] <- soil_profile$soil$LL15[3] * 0.8
  
  
  # Calculate SAT based on porosity and texture
  # note that this is the USA texture triangle - will need to update to the Australian one later
  txt_clss <- texture_class(soil_profile$soil$ParticleSizeClay * 1e-2, soil_profile$soil$ParticleSizeSilt * 1e-2)
  soil_profile$soil$SAT <- calculate_SAT(txt_clss, soil_profile$soil$BD)
  # Ensure DUL is less than SAT
  soil_profile$soil$DUL <- ifelse(soil_profile$soil$DUL > soil_profile$soil$SAT, 
                                  soil_profile$soil$SAT * 0.9, 
                                  soil_profile$soil$DUL)
  # Calculate KS
  B <- (log(1500) - log(33))/(log(soil_profile$soil$DUL) - log(soil_profile$soil$LL15))
  Lambda <- 1/B
  soil_profile$soil$KS <- (1930 * (soil_profile$soil$SAT - soil_profile$soil$DUL)^(3 - Lambda)) * 100
  
  # compute crop-specific LL, KL and XF
  crop_df <- compute_all_crops_params(soil_profile$crops, thickness = soil_profile$soil$Thickness, bd_fine = soil_profile$soil$BD,
                                      intl_sand_propn = soil_profile$soil$ParticleSizeSand/100, 
                                      intl_clay_propn = soil_profile$soil$ParticleSizeClay/100,
                                      ll15 = soil_profile$soil$LL15,
                                      soil_profile$soil$DUL)
  
  common_cols <- intersect(names(soil_profile$soil), names(crop_df))
  soil_profile$soil[common_cols] <- crop_df[common_cols]
  
  
  #### Passing parameters from soilwat
  
  # Get soil parameters
  t2sp <- texture2soilParms(txt_clss)
  
  if(missing(soil.profile)){
    soil_profile$soilwat <- soilwat_parms(Salb = t2sp$Albedo, CN2Bare = t2sp$CN2, 
                                          SWCON = t2sp$SWCON,
                                          Thickness = thcknss)    
  }else{
    soil_profile$soilwat <- soilwat_parms(Salb = t2sp$Albedo, CN2Bare = t2sp$CN2, 
                                          SWCON = t2sp$SWCON,
                                          Thickness = soil_profile$soil$Thickness)    
  }
  
  ### Passing the initial soil water?
  if(new.soil){
    ini.wat <- soil_profile$soil$LL15
    isw <- initialwater_parms(Thickness = soil_profile$soil$Thickness, 
                              InitialValues = ini.wat)
  }else{
    ini.wat <- soil_profile$soil$LL15
    isw <- initialwater_parms(Thickness = thcknss, InitialValues = ini.wat)    
  }
  
  soil_profile$initialwater <- isw
  
  if(find.location.name){
    if(requireNamespace("maps", quietly = TRUE)){
      country <- maps::map.where(x = lon, y = lat)
      if(country == "USA"){
        state <- toupper(maps::map.where(database = "county", x = lon, y = lat)) 
      }else{
        url <- paste0("https://photon.komoot.io/reverse?lon=", lon, "&lat=", lat)
        fgeo <- jsonlite::fromJSON(url)
        state <- fgeo$feature$properties$state
      }
    }else{
      url <- paste0("https://photon.komoot.io/reverse?lon=", lon, "&lat=", lat)
      fgeo <- jsonlite::fromJSON(url)
      state <- fgeo$feature$properties$state
      country <- fgeo$features$properties$country
    }    
  }else{
    state <- ""
    country <- ""
  }
  
  #### Attributes ####
  alist <- list()
  alist$SoilType <- paste("SoilType = ", txt_clss[1])
  alist$State <- state
  alist$Country <- country
  alist$Longitude <- lon
  alist$Latitude <- lat
  alist$DataSource <- paste("Original source is Soil and Landscape Grid of Australia. See: https://esoil.io/TERNLandscapes/Public/Pages/SLGA/index.html",Sys.time())
  alist$Comments <- paste("resolution = 90 m",
                          "- taxonomic classification name =", txt_clss[1],
                          "- drainage class =", NA, 
                          "- elevation =", NA,
                          "- slope =", NA,
                          "- geomdesc =", NA)
  
  soil_profile$metadata <- alist
  
  
  return(soil_profile)  
}


soil_json <- function(template, soil.profile, soil.name){
  # this function is adapted from apsimx package function edit_apsimx_replace_soil_profile.R
  # it set up a new soil profile in json format
  template_j <- template
  
  soil.node0 <- template_j$Children[[1]]$Children
  
  ## First edit: soil 'Physical'
  ## Depth, Thickness, ParticleSizeClay,
  ## BD, AirDry, LL15, DUL, SAT, KS
  wspn <- grepl("Models.Soils.Physical", soil.node0)
  soil.physical.node <- soil.node0[wspn][[1]]
  for(i in c("Depth", "Thickness", "ParticleSizeClay", "ParticleSizeSilt", "ParticleSizeSand", "BD", "AirDry", "LL15", "DUL", "SAT", "KS")){
    ## Format the variable
    if(is.null(soil.profile$soil[[i]])) next
    tmp <- as.vector(soil.profile$soil[[i]], mode = "list")
    ## Replace the variable
    soil.physical.node[[i]] <- tmp 
  }
  ## Preliminary setup of soil
  if(length(soil.profile$crops) > length(soil.physical.node$Children)){
    for(i in seq_along(soil.profile$crops)){
      soil.physical.node$Children[[i]] <- soil.physical.node$Children[[1]]
      soil.physical.node$Children[[i]]$Name <- paste0(soil.profile$crops[i], "Soil")
    }
  }
  ## Crop parameters
  for(i in 1:length(soil.profile$crops)){
    crop.vars.names <- paste0(soil.profile$crops[i], c(".KL", ".LL", ".XF"))
    for(j in crop.vars.names){
      tmp <- as.vector(soil.profile$soil[[j]], mode = "list")
      strpcrop <- strsplit(j, ".", fixed = TRUE)[[1]][2]
      soil.physical.node$Children[[i]][[strpcrop]] <- tmp
    }
  }
  soil.node0[wspn][[1]] <- soil.physical.node
  template_j$Children[[1]]$Children <- soil.node0
  
  ## Next edit the soil organic component
  wsomn <- grepl("Models.Soils.Organic", soil.node0)
  soil.om.node <- soil.node0[wsomn][[1]]
  
  for(i in c("Depth", "Thickness", "Carbon", "SoilCNRatio", "FBiom", "FInert", "FOM")){
    ## Format the variable
    tmp <- as.vector(soil.profile$soil[[i]], mode = "list")
    ## Replace the variable
    soil.om.node[[i]] <- tmp 
  }
  soil.node0[wsomn][[1]] <- soil.om.node
  template_j$Children[[1]]$Children <- soil.node0
  
  ## Next edit the Chemical component
  wschn <- grepl("Models.Soils.Chemical", soil.node0)
  soil.chemical.node <- soil.node0[wschn][[1]]
  
  for(i in c("Depth", "Thickness", "PH", "CEC")){
    ## Format the variable
    tmp <- as.vector(soil.profile$soil[[i]], mode = "list")
    ## Replace the variable
    soil.chemical.node[[i]] <- tmp 
  }
  soil.node0[wschn][[1]] <- soil.chemical.node
  template_j$Children[[1]]$Children <- soil.node0
  
  ## Next edit the Solute component, this is not present in older versions of APSIM
  wssoln <- grepl("Models.Soils.Solute", soil.node0)
  if(sum(wssoln) > 0){
    soil.solute.nodes <- soil.node0[wssoln]
    for(i in seq_along(soil.solute.nodes)){
      soil.solute.node <- soil.solute.nodes[[i]]
      
      vname <- soil.solute.node$Name
      if (!is.null(vname)){
        vname.in.soil <- grep(vname, names(soil.profile$soil), value = TRUE)
        if(length(vname.in.soil) > 0){
          tmp <- as.vector(soil.profile$soil[[vname.in.soil]], mode = "list") ## This now works because these are initial values?
          soil.solute.node$InitialValues <- tmp  
          soil.solute.node$Thickness <- soil.profile$soil$Thickness
          soil.node0[wssoln][[i]] <- soil.solute.node
        }
      }
      
    }
    template_j$Children[[1]]$Children <- soil.node0    
  }
  
  ## Edit metadata
  if(!is.null(soil.profile$metadata)){
    soil.node.names <- names(template_j$Children[[1]])
    skp.nms <- c("$type", "Name", "Children", "IncludeInDocumentation", "Enabled", "ReadOnly")
    for(i in soil.node.names){
      if(i %in% skp.nms) next
      if(i %in% names(soil.profile$metadata)){
        template_j$Children[[1]][[i]] <- soil.profile$metadata[[i]]
      }else{
        if(i %in% c("RecordNumber","ApsoilNumber")){
          template_j$Children[[1]][[i]] <- 0 
        }else{
          template_j$Children[[1]][[i]] <- NA  
        }
      }
    }
  }
  
  ## Use metadata to edit the name of the soil
  
  template_j$Children[[1]]$Name <- soil.name
  
  ## Edit soilWat if present
  if(inherits(soil.profile$soilwat, "soilwat_parms")){
    
    wsswn <- grepl("Models.WaterModel", soil.node0) 
    soil.soilwat.node <- soil.node0[wsswn][[1]]
    for(i in names(soil.profile$soilwat)){
      if(any(is.na(soil.profile$soilwat[[i]]))){
        next
      }else{
        soil.soilwat.node[[i]] <- soil.profile$soilwat[[i]]  
      } 
    }
    soil.node0[wsswn][[1]] <- soil.soilwat.node
    template_j$Children[[1]]$Children <- soil.node0
  }
  
  ## Edit InitialWater if present
  if(inherits(soil.profile$initialwater, "initialwater_parms")){
    
    wsiswn <- grepl("Models.Soils.Water", soil.node0) 
    soil.initialwater.node <- soil.node0[wsiswn][[1]]
    for(i in names(soil.profile$initialwater)){
      if(any(is.na(soil.profile$initialwater[[i]]))){
        next
      }else{
        soil.initialwater.node[[i]] <- soil.profile$initialwater[[i]]  
      } 
    }
    soil.node0[wsiswn][[1]] <- soil.initialwater.node
    template_j$Children[[1]]$Children <- soil.node0
  }
  
  ## Edit Solutes if present
  if(inherits(soil.profile$solutes, "solutes_parms")){
    
    wssoln <- grepl("Models.Soils.Solute", soil.node0) 
    soil.solute.nodes <- soil.node0[wssoln]
    for(i in seq_along(soil.solute.nodes)){
      soil.solute.node <- soil.solute.nodes[[i]]
      vname <- soil.solute.node$Name
      vname.in.solutes <- grep(vname, soil.profile$solutes$Solutes, value = TRUE)
      if(length(vname.in.solutes) > 0){
        for(j in seq_along(soil.profile$solutes)){
          if(any(is.na(soil.profile$solutes[[j]]))){
            next
          }else{
            sspntc <- names(soil.profile$solutes)[[j]] ## Soil Solute Parameter Name to Change
            ## Is the parameter above present in soil.solute.node?
            soil.solute.node.names <- names(soil.solute.node)
            if(!sspntc %in% soil.solute.node.names)
              stop("soil solute parameter to change is not in soil.solute.node")
            soil.solute.node[sspntc] <- soil.profile$solutes[[j]]    
          } 
        }        
      }
      soil.solute.nodes[[i]] <- soil.solute.node
    }
    soil.node0[wssoln] <- soil.solute.nodes
    template_j$Children[[1]]$Children <- soil.node0
  }
  return(template_j$Children[[1]])
  
}

# example script --------------------

json_file <- "soil_template.apsimx"  # Path to the JSON file
template_json <- read_json(json_file)
skeleton_json <- template_json

crops_df <- data.frame(crops = c("Wheat", "Sorghum", "Chickpea", "Mungbean",
                                 "Canola", "Barley", "SCRUM",
                                 "Oats", "Lucerne"))

template_soil <- apsimx_soil_profile(Thickness = c(150, 150, 150, 150, 200, 200, 200, 200, 300, 300),
                                     crops = c("Wheat", "Sorghum", "Chickpea", "Mungbean",
                                               "Canola", "Barley", "SCRUM",
                                               "Oats", "Lucerne"))

# the location used show a reduction of xf factor down the depth due to bulk density and sand/clay content

# this obtains the soil profile conforming to the depth intervals used in the World modellers database, 
# where kl for the last depth (i.e., 170-200 cm) is reduced to 0 due to the high clay content 
# (rooting depths decreased in the function by 30 cm for soils with clay% > 40 throughout the profile)
WM_depth <- get_slga_soil_profile(c(150.397, -28.127), soil.profile = template_soil)

# this obtains the soil profile using native slga depths 
slga_depth <- get_slga_soil_profile(c(150.397, -28.127), xargs = crops_df)

WM_depth_j <- soil_json(template_json, WM_depth, "SLGA_WM_depth")
slga_depth_j <- soil_json(template_json, slga_depth, "SLGA_default_depth")

skeleton_json$Children[[1]] <- WM_depth_j
skeleton_json$Children[[2]] <- slga_depth_j

write_json(skeleton_json, "slga_test.apsimx", 
           pretty = TRUE, digits = NA, auto_unbox = TRUE, , null = "null",
           na = "null")
