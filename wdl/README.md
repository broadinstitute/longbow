# MAS-seq / Longbow Workflows

This folder contains the workflow for processing MAS-seq data with
Longbow and all other required tools to produce a transcript count
matrix for further analysis.

These [WDL](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md) files have several numbered tasks
which are executed in approximate numbered order when submitted to the cloud via [Cromwell](https://cromwell.readthedocs.io/en/stable/) or [Terra](https://terra.bio/)
(the order is not exactly numerical because cromwell performs call tree analysis to 
parallelize tasks that can be run concurrently). 

## Alpha Warning
The workflows defined here should be considered an alpha release.
They produce valid data, however they are not optimized and are
therefore more expensive than they need to be.

## Docker Images
All docker images used in these WDL files are stored in the 
[Google Container Registry](https://cloud.google.com/container-registry) for the DSP Methods Long Reads Methods and 
Applications group at The Broad Institute of MIT and Harvard.

