############### MSOM code for WA_ExpertReview_vs_AIReview in WA study area:

############### Initial Code to create edited model text with added BPV Mu_Psi and Mu_p 

## Load packages

library(camtrapR)
library(purrr)
library(DT)
library(knitr)
library(ggplot2)
library(sf)
library(rlist)
library(rjags)
library(coda)
library(gridExtra)
library(dplyr)


### working directory | need to adjust to study

setwd("C:/Users/Home/Documents/2024_5_Fall_PhD_Work/AI_ManuscriptWork/Lynx2016-17_MSOM_V2")

getwd()

### call camop, record table, site covs

camtraps <- read.csv("2016_2017_OK_KET_Master_CameraInfo_cleaned.csv")

sitecovs <-  read.csv("Lynx_2016_2017_AI_Covariates_Cleaned.csv")

##### only table changing between tests:

recordTable <- read.csv("2016_2017_OK_KET_60minDelay_recordtable_V4_SpeciesLimited.csv")


#########
# Covariate preprocessing | all categorical or non-numeric covariates must be made into a factor or numeric value

sitecovs$period <- as.factor(sitecovs$period) 
sitecovs$trailtype <- as.factor(sitecovs$trailtype) 
sitecovs$Angle <- as.factor(sitecovs$Angle) 
sitecovs$Kettles_or_Okanogan <- as.factor(sitecovs$Kettles_or_Okanogan) 

######## Scale all covariates of interest
cols_to_scale <- c("Duration", "JulianDate_Deployment", "ghm", "d2hs", "elevation", "slope", "aspect", 
                   "canopy", "evi", "precipitation", "temperature")
sitecovs <- sitecovs %>%
  mutate(across(all_of(cols_to_scale), ~ as.numeric(scale(.x, center = TRUE, scale = TRUE))))



############ running a no data NA camop code or other function breaks in deployed cameras.| camera operation matrix
 camop_no_problem <- cameraOperation(CTtable      = camtraps,
                                    stationCol   = "Station",
                                    setupCol     = "Setup_date",
                                    retrievalCol = "Retrieval_date.Function.Date",
                                    hasProblems  = FALSE,
                                    dateFormat   = "%m/%d/%Y"
)

####Multi-species combined detection histories | Full cameras; Full Species; 10-Occ Length (adjust as appropriate)

################ creating the detection history package from record table. Here I need to specify occasion length, specify record table structure, and if you run detection history as start date of station or start date of study.

 DetHist_list <- lapply(unique(recordTable$Species), FUN = function(x) {
   detectionHistory(
     recordTable         = recordTable,
     camOp                = camop_no_problem,
     stationCol           = "Station",
     speciesCol           = "Species",
     recordDateTimeCol    = "DateTimeOriginal",
     recordDateTimeFormat = "%m/%d/%Y %H:%M",
     species              = x,     # this gets modifies by lapply
     occasionLength       = 10,
     day1                 = "station",
     maxNumberDays        =  150,
     datesAsOccasionNames = FALSE,
     includeEffort        = TRUE,
     scaleEffort          = FALSE,
     timeZone             = "UTC"
   )}
 )

 # assign species names to the list items

 names(DetHist_list) <- unique(recordTable$Species)

##Get the detection history of each species and put into a new list (thereby removing the effort matrix).

 ylist <- lapply(DetHist_list, FUN = function(x) x$detection_history)

####### Create Data List  bundle the necessary data for communityModel

 data_list <- list(ylist    = ylist,
                   siteCovs = sitecovs,
                   obsCovs  = list(effort = DetHist_list[[1]]$effort))



####### Build model structure and final data package:

# temp model text used here as we will extract the text and add in the model diagnostics for the full model run:

  modelfile <- tempfile(fileext = ".txt")

  mod.jags <- communityModel(
    data_list,
    occuCovs = list(ranef = c("d2hs", "elevation", "slope", "canopy", "evi", "precipitation")),
    detCovs = list(ranef = c("Duration", "JulianDate_Deployment")),
    detCovsObservation = list(fixed = "effort"),
    intercepts = list(det = "ranef", occu = "ranef"),
    modelFile = modelfile
  )

############# follow temp model save location and extract model text code.

############# add corrected calculations for BPV and Mu_PSI Mu_P

############## See new model text 'Lynx2016_2017_MSOM_T6_modeltext_ParametersAdded.txt'


####################################################################################################################################################################################################################################################


############## full MSOM code run with edited model text on the Kamiak High Preformance Computer at WSU 

## Load packages:

library(camtrapR)
library(purrr)
library(DT)
library(knitr)
library(ggplot2)
library(sf)
library(rlist)
library(rjags)
library(coda)
library(gridExtra)
library(dplyr)


# Static inputs
camtraps <- read.csv("2016_2017_OK_KET_Master_CameraInfo_cleaned.csv")
sitecovs <- read.csv("Lynx_2016_2017_AI_Covariates_Cleaned.csv")

# Factorize categorical covariates
sitecovs$period <- as.factor(sitecovs$period)
sitecovs$trailtype <- as.factor(sitecovs$trailtype)
sitecovs$Angle <- as.factor(sitecovs$Angle)
sitecovs$Kettles_or_Okanogan <- as.factor(sitecovs$Kettles_or_Okanogan)

# Scale numerical covariates
cols_to_scale <- c("Duration", "JulianDate_Deployment", "ghm", "d2hs", "elevation", "slope", "aspect", 
                   "canopy", "evi", "precipitation", "temperature")
sitecovs <- sitecovs %>%
  mutate(across(all_of(cols_to_scale), ~ as.numeric(scale(.x, center = TRUE, scale = TRUE))))

# Create camera operation matrix
camop_no_problem <- cameraOperation(
  CTtable = camtraps,
  stationCol = "Station",
  setupCol = "Setup_date",
  retrievalCol = "Retrieval_date.Function.Date",
  hasProblems = FALSE,
  dateFormat = "%m/%d/%Y"
)

# Record table file names
record_tables <- c("2016_2017_OK_KET_60minDelay_recordtable_V4_SpeciesFull.csv",
  "2016_2017_OK_KET_60minDelay_recordtable_V4_SpeciesLimited.csv",
  "lynx_2016_2017_FullSpeciesCleanedPairs_Class50_recordTable_V4_60min.csv",
  "lynx_2016_2017_FullSpeciesCleanedPairs_Class85_recordTable_V4_60min.csv",
  "lynx_2016_2017_FullSpeciesCleanedPairs_Class95_recordTable_V4_60min.csv",
  "lynx_2016_2017_FullSpeciesRaw_Class50_recordTable_V4_60min.csv",
  "lynx_2016_2017_FullSpeciesRaw_Class85_recordTable_V4_60min.csv",
  "lynx_2016_2017_FullSpeciesRaw_Class95_recordTable_V4_60min.csv",
  "lynx_2016_2017_LimitedSpeciesCleanedPairs_Class50_recordTable_V4_60min.csv",
  "lynx_2016_2017_LimitedSpeciesCleanedPairs_Class85_recordTable_V4_60min.csv",
  "lynx_2016_2017_LimitedSpeciesCleanedPairs_Class95_recordTable_V4_60min.csv",
  "lynx_2016_2017_LimitedSpeciesRaw_Class50_recordTable_V4_60min.csv",
  "lynx_2016_2017_LimitedSpeciesRaw_Class85_recordTable_V4_60min.csv",
  "lynx_2016_2017_LimitedSpeciesRaw_Class95_recordTable_V4_60min.csv" 
)

# Flexible datetime format detection
detect_datetime_format <- function(datetime_vec) {
  x <- head(na.omit(datetime_vec), 10)
  formats <- c(
    "%Y-%m-%d %H:%M:%S",   # ISO format
    "%m/%d/%Y %H:%M:%S",   # US format, with seconds
    "%m/%d/%Y %H:%M"       # US format, no seconds
  )
  for (fmt in formats) {
    parsed <- suppressWarnings(as.POSIXct(x, format = fmt, tz = "UTC"))
    if (all(!is.na(parsed))) return(fmt)
  }
  stop("Unrecognized datetime format in 'DateTimeOriginal'.")
}

# Main loop through record tables
for (rec_file in record_tables) {
  message("Running MSOM for: ", rec_file)
  
  prefix <- tools::file_path_sans_ext(rec_file)
  recordTable <- read.csv(rec_file)

  # Detect appropriate datetime format
  datetime_format <- detect_datetime_format(recordTable$DateTimeOriginal)

  # Build detection histories
  DetHist_list <- lapply(unique(recordTable$Species), FUN = function(x) {
    detectionHistory(
      recordTable = recordTable,
      camOp = camop_no_problem,
      stationCol = "Station",
      speciesCol = "Species",
      recordDateTimeCol = "DateTimeOriginal",
      recordDateTimeFormat = datetime_format,
      species = x,
      occasionLength = 10,
      day1 = "station",
      maxNumberDays = 150,
      datesAsOccasionNames = FALSE,
      includeEffort = TRUE,
      scaleEffort = FALSE,
      timeZone = "UTC"
    )
  })
  names(DetHist_list) <- unique(recordTable$Species)
  ylist <- lapply(DetHist_list, function(x) x$detection_history)

  # Create data bundle
  data_list <- list(
    ylist = ylist,
    siteCovs = sitecovs,
    obsCovs = list(effort = DetHist_list[[1]]$effort)
  )

  # Build and fit model
  modelfile <- tempfile(fileext = ".txt")
  mod.jags <- communityModel(
    data_list,
    occuCovs = list(ranef = c("d2hs", "elevation", "slope", "canopy", "evi", "precipitation")),
    detCovs = list(ranef = c("Duration", "JulianDate_Deployment")),
    detCovsObservation = list(fixed = "effort"),
    intercepts = list(det = "ranef", occu = "ranef"),
    modelFile = modelfile
  )

############ Import in the edited text file of the model structure (BPV Psi and P recording)

modelFile_edited = 'Lynx2016_2017_MSOM_T6_modeltext_ParametersAdded.txt'

############### add the additional parameters to watch in model loop

params = c("mu.psi", "mu.p", "Bpvalue", "Bpvalue_species", mod.jags@params)

## full run numbers
n.iter = 40000
n.burnin = 20000
thin = 50
chains = 5

########

mod <- rjags::jags.model(file =  modelFile_edited, 
                               data =  mod.jags@data, 
                               inits =  mod.jags@inits_fun(),
                               n.chain=chains, 
                               n.adapt=0,
                               quiet = TRUE)
      
      out <- rjags::coda.samples(model = mod,
                                variable.names =  params, 
                                n.iter	= n.iter, 
                                thin = thin)


   out_mcmclist <- coda::mcmc.list(out)

   fit.jags <- window(out_mcmclist, 
                     start=n.burnin+1, 
                     end = n.iter)


  # Save model and fit
  list.save(mod.jags, paste0(prefix, "_mod2.rds"))
  list.save(fit.jags, paste0(prefix, "_fit2.rds"))

}





