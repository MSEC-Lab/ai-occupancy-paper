# Overview

This repository contains supplementary code for the paper "Review of camera trap images by experts and artificial intelligence yield similar multi-species occupancy models".

## Image processing

The image processing described in the paper follows [this notebook](https://github.com/agentmorris/MegaDetector/blob/main/notebooks/manage_local_batch.py), which performs the following steps:

* Runs [MegaDetector](https://github.com/agentmorris/MegaDetector/) on images
* Runs [SpeciesNet](https://github.com/google/cameratrapai) on the animals detected by MegaDetector
* Applies the taxonomic mapping described in the paper, via the [restrict_to_taxa_list](https://megadetector.readthedocs.io/en/latest/postprocessing.html#megadetector.postprocessing.classification_postprocessing.restrict_to_taxa_list) function
* Applies the prediction smoothing described in the paper, via the [smooth_classification_results_image_level](https://megadetector.readthedocs.io/en/latest/postprocessing.html#megadetector.postprocessing.classification_postprocessing.smooth_classification_results_image_level) and [smooth_classification_results_sequence_level](https://megadetector.readthedocs.io/en/latest/postprocessing.html#megadetector.postprocessing.classification_postprocessing.smooth_classification_results_image_level) functions

