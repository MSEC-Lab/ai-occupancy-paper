############### MSOM code for GT_ExpertReview_vs_AIReview in GT study area:

############### Initial Code to created edited model text with added BPV Mu_Psi and Mu_p 

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

setwd("~/1_AI_MSOMs/Guatamala_MSOM_T6/Mazama_Added")

getwd()


################## Read all needed files --------------

# Load static inputs
camtraps <- read.csv("CamOp_DryOnly_cleaned.csv")
sitecovs <- read.csv("covariates_final_scaled_dry_clean.csv")

##### only table changing between tests:

recordTable <- read.csv("RT_raw95_60min_ltd_DryOnly_MSp.csv")

# Factorize categorical covariates
# Covariate preprocessing | all categorical or non-numeric covariates must be made into a factor or numeric value

sitecovs$feature <- as.factor(sitecovs$feature)
sitecovs$Station <- as.factor(sitecovs$Station)

### note covariate file already scaled so no additional scaling needed

######## Create a camera operation matrix using Camop file

camop_no_problem <- cameraOperation(
  CTtable = camtraps,
  stationCol = "Station",
  setupCol = "Setup_date",
  retrievalCol = "Retrieval_date.Function.Date",
  hasProblems = FALSE,
  dateFormat = "%m/%d/%Y %H:%M"
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
     occasionLength       = 3,
     day1                 = "station",
     maxNumberDays        =  60,
     datesAsOccasionNames = FALSE,
     includeEffort        = TRUE,
     scaleEffort          = FALSE,
     timeZone             = "America/Guatemala"
   )}
 )

 # assign species names to the list items

 names(DetHist_list) <- unique(recordTable$Species)

##Get the detection history of each species and put into a new list (thereby removing the effort matrix).

 ylist <- lapply(DetHist_list, FUN = function(x) x$detection_history)

#### Create Data List  bundle the necessary data for communityModel

 data_list <- list(ylist    = ylist,
                   siteCovs = sitecovs,
                   obsCovs  = list(effort = DetHist_list[[1]]$effort))

####### Build model structure and final data package:

# temp model text used here as we will extract the text and add in the model diagnostics for the full model run:

  modelfile <- tempfile(fileext = ".txt")
 
mod.jags <- communityModel(
        data_list,
        occuCovs = list(ranef = c("d2hs", "d2hs_squared", 
                                  "elevation", "slope", "d2w",
                                  "canopy", "evi", "precipitation")),
        detCovs = list(ranef = c("feature","duration", "jDate")),
        detCovsObservation = list(fixed = "effort"),
        intercepts = list(det = "ranef", occu = "ranef"),
        modelFile = modelfile) 
      
   

############# follow temp model save location and extract model text code.

############# add corrected calculations for BPV and Mu_PSI Mu_P

############## See new model text 'Guatemala_MSOM_T6_modeltext_ParametersAdded.txt'


####################################################################################################################################################################################################################################################



############## full MSOM code run with edited model text on the Kamiak High Preformance Commuter at WSU 

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
camtraps <- read.csv("CamOp_DryOnly_cleaned.csv")
sitecovs <- read.csv("covariates_final_scaled_dry_clean.csv")


# Factorize categorical covariates
sitecovs$feature <- as.factor(sitecovs$feature)
sitecovs$Station <- as.factor(sitecovs$Station)

# Camera operation matrix (no problem periods used)
camop_no_problem <- cameraOperation(
  CTtable = camtraps,
  stationCol = "Station",
  setupCol = "Setup_date",
  retrievalCol = "Retrieval_date.Function.Date",
  hasProblems = FALSE,
  dateFormat = "%m/%d/%Y %H:%M"
)

# Input record tables
record_tables <- c(
  "RT_CleanedPairs50_60min_full_DryOnly_MSp.csv",
  "RT_CleanedPairs50_60min_ltd_DryOnly_MSp.csv",
  "RT_CleanedPairs85_60min_full_DryOnly_MSp.csv",
  "RT_CleanedPairs85_60min_ltd_DryOnly_MSp.csv",
  "RT_CleanedPairs95_60min_full_DryOnly_MSp.csv",
  "RT_CleanedPairs95_60min_ltd_DryOnly_MSp.csv",
  "RT_expertReview_full_DryOnly_MSp.csv",
  "RT_expertReview_ltd_DryOnly_MSp.csv",
  "RT_raw50_60min_full_DryOnly_MSp.csv",
  "RT_raw50_60min_ltd_DryOnly_MSp.csv",
  "RT_raw85_60min_full_DryOnly_MSp.csv",
  "RT_raw85_60min_ltd_DryOnly_MSp.csv",
  "RT_raw95_60min_full_DryOnly_MSp.csv",
  "RT_raw95_60min_ltd_DryOnly_MSp.csv"
)

# Main loop through record tables
for (rec_file in record_tables) {
  message("Running MSOM for: ", rec_file)
  
  prefix <- tools::file_path_sans_ext(rec_file)
  recordTable <- read.csv(rec_file)
 
 # Build detection histories

DetHist_list <- lapply(unique(recordTable$Species), FUN = function(x) {
   detectionHistory(
     recordTable         = recordTable,
     camOp                = camop_no_problem,
     stationCol           = "Station",
     speciesCol           = "Species",
     recordDateTimeCol    = "DateTimeOriginal",
     recordDateTimeFormat = "%m/%d/%Y %H:%M",
     species              = x,     # this gets modifies by lapply
     occasionLength       = 3,
     day1                 = "station",
     maxNumberDays        =  60,
     datesAsOccasionNames = FALSE,
     includeEffort        = TRUE,
     scaleEffort          = FALSE,
     timeZone             = "America/Guatemala"
   )}
 )


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
        occuCovs = list(ranef = c("d2hs", "d2hs_squared", 
                                  "elevation", "slope", "d2w",
                                  "canopy", "evi", "precipitation")),
        detCovs = list(ranef = c("feature","duration", "jDate")),
        detCovsObservation = list(fixed = "effort"),
        intercepts = list(det = "ranef", occu = "ranef"),
        modelFile = modelfile) 

 ############ import in a saved text file of the model structure from the tmp model text file on main computer which added BPV Psi and P recording
  
  modelFile_edited = 'Guatemala_MSOM_T6_modeltext_ParametersAdded.txt'
  
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
  list.save(mod.jags, paste0(prefix, "_mod.rds"))
  list.save(fit.jags, paste0(prefix, "_fit.rds"))
  
}


