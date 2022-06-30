# Overview

This script is used to obtain the presence proportion of study species and the occurrence proportion of landcover types. This script can only be run after first successfully running the preprocessing code. There is no need to move the processed data from its location within the *1.preprocessing* folder structure.

Note: this code is fast and so can be easily run on a local machine, provided you have the processed data required. 

## Inputs 
1. Processed eBird checklist data for all species.
2. At least one processed covariate prediction file. 

## Outputs 
1. Presence proportions for hunted or non-hunted study species: `sighting_props.csv`
2. Occurrence proportions for all landcover types: `landcover_props.csv`
