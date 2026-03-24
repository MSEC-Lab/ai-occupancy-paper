# Overview

This repository contains supplementary code for the paper "Rapid ecological inference for megadiverse communities: Identification of camera trap images by artificial intelligence and human experts produce similar multi-species occupancy models".

## Image processing

The image processing described in the paper follows [this notebook](https://github.com/agentmorris/MegaDetector/blob/main/notebooks/manage_local_batch.py), which performs the following steps:

* Runs [MegaDetector](https://github.com/agentmorris/MegaDetector/) on images
* Runs [SpeciesNet](https://github.com/google/cameratrapai) on the animals detected by MegaDetector
* Applies the taxonomic mapping described in the paper, via the [restrict_to_taxa_list](https://megadetector.readthedocs.io/en/latest/postprocessing.html#megadetector.postprocessing.classification_postprocessing.restrict_to_taxa_list) function
* Applies the prediction smoothing described in the paper, via the [smooth_classification_results_image_level](https://megadetector.readthedocs.io/en/latest/postprocessing.html#megadetector.postprocessing.classification_postprocessing.smooth_classification_results_image_level) and [smooth_classification_results_sequence_level](https://megadetector.readthedocs.io/en/latest/postprocessing.html#megadetector.postprocessing.classification_postprocessing.smooth_classification_results_image_level) functions

## Data analysis

* [Guatemala_MSOM.R](Guatemala_MSOM.R) runs the Bayesian multi-species occupancy model with covariates for the Guatemala study area.
* [Montana_MSOM.R](Montana_MSOM.R) runs the Bayesian multi-species occupancy model with covariates for the Montana study area.
* [Washington_MSOM.R](Washington_MSOM.R) runs the Bayesian multi-species occupancy model with covariates for the Washington study area.
* [Spatial_projection.R](Spatial_projection.R) projects occupancy model across each study area for AI and expert MSOMs, and calculates "difference" maps.

## Data files
For Washington study site, there are a number of files that are referred to in the Washington_MSOM code:
* 2016_2017_OK_KET_Master_CameraInfo_cleaned.csv - Camera operation table (Exludes Latitude/Longitude)
* Lynx_2016_2017_AI_Covariates_Cleaned.csv - Covariate file 
* 2016_2017_OK_KET_60minDelay_recordtable_V4_SpeciesFull.csv  - Record table of species detections for expert reviewed, full dataset
* 2016_2017_OK_KET_60minDelay_recordtable_V4_SpeciesLimited.csv  - Record table of species detections for expert reviewed, limited species dataset
* lynx_2016_2017_FullSpeciesCleanedPairs_Class50_recordTable_V4_60min.csv  - Record table of species detections for AI reviewed, full dataset with 50% confidence cutoff and cleaned sequences
* lynx_2016_2017_FullSpeciesCleanedPairs_Class85_recordTable_V4_60min.csv  - Record table of species detections for AI reviewed, full dataset with 85% confidence cutoff and cleaned sequences
* lynx_2016_2017_FullSpeciesCleanedPairs_Class95_recordTable_V4_60min.csv  - Record table of species detections for AI reviewed, full dataset with 95% confidence cutoff and cleaned sequences
* lynx_2016_2017_FullSpeciesRaw_Class50_recordTable_V4_60min.csv  - Record table of species detections for AI reviewed, full dataset with 50% confidence cutoff and uncleaned sequences
* lynx_2016_2017_FullSpeciesRaw_Class85_recordTable_V4_60min.csv  - Record table of species detections for AI reviewed, full dataset with 85% confidence cutoff and uncleaned sequences
* lynx_2016_2017_FullSpeciesRaw_Class95_recordTable_V4_60min.csv  - Record table of species detections for AI reviewed, full dataset with 95% confidence cutoff and uncleaned sequences
* lynx_2016_2017_LimitedSpeciesCleanedPairs_Class50_recordTable_V4_60min.csv  - Record table of species detections for AI reviewed, limited dataset with 50% confidence cutoff and cleaned sequences
* lynx_2016_2017_LimitedSpeciesCleanedPairs_Class85_recordTable_V4_60min.csv  - Record table of species detections for AI reviewed, limited dataset with 85% confidence cutoff and cleaned sequences
* lynx_2016_2017_LimitedSpeciesCleanedPairs_Class95_recordTable_V4_60min.csv  - Record table of species detections for AI reviewed, limited dataset with 95% confidence cutoff and cleaned sequences
* lynx_2016_2017_LimitedSpeciesRaw_Class50_recordTable_V4_60min  - Record table of species detections for AI reviewed, limited dataset with 50% confidence cutoff and uncleaned sequences
* lynx_2016_2017_LimitedSpeciesRaw_Class85_recordTable_V4_60min  - Record table of species detections for AI reviewed, limited dataset with 85% confidence cutoff and uncleaned sequences
* lynx_2016_2017_LimitedSpeciesRaw_Class95_recordTable_V4_60min  - Record table of species detections for AI reviewed, limited dataset with 95% confidence cutoff and uncleaned sequences
## Data availability

A subset of the images used for this paper are available as the "[WSU Lynx](https://lila.science/datasets/wsu-lynx/)" dataset on [LILA BC](https://lila.science).
