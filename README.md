
**Code to reproduce:**

McFadden et al. 2019. Disentangling the functional trait correlates of spatial aggregation in tropical forest trees. Ecology

<br/>

**Note:**

This code and the analyses it performs largely follow the approach of Shen et al. (2013, see Appendix B and code citation below), except this study used a negative exponential distribution for the shape parameter of the Matern covariance function (0.5) whereas in Shen et al. 2013 this was a free parameter that was fit from a range of values. In addition, the proportion of variance explained by habitat and dispersal within species was not estimated. 


*References and Code Cited*

Shen, G., F. He, R. Waagepetersen, I.-F. Sun, Z. Hao, Z.-S. Chen, and M. Yu. 2013. Quantifying effects of habitat heterogeneity and other clustering processes on spatial distributions of tree species. Ecology 94:2436-2443.

Shen, G., F. He, R. Waagepetersen, I.-F. Sun, Z. Hao, Z.-S. Chen, and M. Yu. 2016, August 10. Supplement 1. R functions for estimating parameters of the Cox process, model selections, and calculating percentage of variances explained by habitat heterogeneity and other clustering processes.
DOI: https://doi.org/10.6084/m9.figshare.3557760.v1


<br/>

**Generalized vignette using the functions contained in this repository to reproduce McFadden et al. 2019**

```r

## NOTE:
## This script uses dummy names for stem and habitat data
## Change the names to your stem data to use the function

## From:
## McFadden et al. 2019. Disentangling the functional trait correlates
## of spatial aggregation in tropical forest trees. Ecology 


### Setup #################=
# Need to make sure both these functions are in the working directory:
# 'sppm_decision_tree.r'
# 'gof_test.R' 

rm(list=ls())

# Source the functions
source("gof_test.R")
source("sppm_decision_tree.r")

# Load stem data
load("cen.Rdata") # Change this to your actual tree stem data file

# Load list of habitat attribute image objects of class spatat::im
load("habList.Rdata") # Change this to actual im object list

# To choose a few low abundance species to run first to test the function remove the #s below 
# abund <- names(which(table(cen$sp) > 70 & table(cen$sp) < 75))
# cen <- cen[cen$sp %in% abund,]
# length(unique(cen$sp)) # How many species were chosen?

# Run the decision tree function
res <- sppm_decision_tree(dataframe=cen,           # Stem coordinate dataset
                          species_codes='sp',      # Name of species code column 
                          x_coords='gx',           # X coordinate column name
                          y_coords='gy',           # Y coordinate column name 
                          hab_vars=habList,        # Habitat attribute im object list name
                          nu=0.5,                  # Preset to exponential decay, can change this if desired 
                          abund_cutoff=70,         # Preset to 70, can change this if desired 
                          plotdim=c(0,500,0,500))  # Plot dimensions in meters, your plot dimensions may differ 

# Print the output 
res

```
