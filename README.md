# Overview

This repository contains supplementary code for the paper "Review of camera trap images by experts and artificial intelligence yield similar multi-species occupancy models".

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
