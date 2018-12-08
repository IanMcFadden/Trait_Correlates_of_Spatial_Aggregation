
## Ian McFadden
## imcfadden@ucla.edu

################################################################################
#
# sppm_decision_tree
#
# Implements a decision tree using spatial point process
# models (SPPMs) to:
#
# 1.
# Categorize species based on spatial aggregation into:
#    Category 1: 'Complete spatial randomness' (CSR)
#    Category 2: 'Habitat only' (H)
#    Category 3: 'Dispersal only' (D)
#    Category 4: 'Habitat and dispersal' (HD)
#
# 2.
# When appropriate estimates model parameters related to habitat 
# association and / or dispersal limitation based on species category
# using SPPMs incorporating one or both processes:
#
#    Category 1 (CSR): No parameters estimated
#
#    Category 2 (H):   Heterogenious poisson process to estimete coefficients
#                      of association only  
#
#    Category 3 (D):   Thomas cluster process to estimate cluster size (alpha)
#                      and clustering intensity (sigma2) only 
#
#    Category 4 (HD):  log-Gaussian Cox process to estimate coefficients
#                      of association, cluster size (alpha) and clustering 
#                      intensity (sigma2)
#
# ARGUMENTS
# dataframe       Spatial dataframe to be used 
# species_codes   Name of species code column in dataframe
# x_coords        Name of x coordinate column in dataframe
# y_coords        Name of y coordinate column in dataframe
# hab_vars        List of habitat images in 'im' format
# nu              Preset kernel decay shape parameter 
#                   0.5 (default): exponential
#                   Inf (infinity): Gaussian      
# abund_cutoff    Species stem abundance cutoff, 70 is default 
# plotdim         Dimensions of the study plot as c(xmin, xmax, ymin, ymax) 
#
# MODEL OUTPUTS
# *Output will differ by category (C1-4)
# Coefficients of association with habitat variables 
# CI95lo          Min of 95% confidence interval 
# CI95hi          Max of 95% confidence interval 
# sigma2          Clustering intensity, the mean number of points per cluster
# alpha           Cluster size, in meters 
#
################################################################################

sppm_decision_tree <- function(dataframe=NULL,
                               species_codes=NULL,
                               x_coords=NULL,
                               y_coords=NULL,
                               hab_vars=NULL,
                               nu=0.5,
                               abund_cutoff=70,
                               plotdim=c(0,500,0,500)){

     # Require packages 
     packages <- c("spatstat", "raster", "geostatsp", "ncf", "maptools",
                   "spatial", "RandomFieldsUtils", "RandomFields")
     lapply(packages, require, character.only=T)
     
     # Source the goodness of fit function 
     source("gof_test.r")
     
     # Remove species below abundance cutoff
     abund <- names(which(table(dataframe[,species_codes]) >= abund_cutoff))
     dataframe <- dataframe[dataframe[,species_codes] %in% abund,]
     
     # Make vector of species to loop over
     spp <- unique(dataframe[,species_codes])
     
     # Make vector of habitat variable names 
     habNames <- names(hab_vars)
     
     # Create column names
     CIlo <- paste0(names(hab_vars), "CI95.lo")
     CIhi <- paste0(names(hab_vars), "CI95.hi")
     colNames <- c("sp", "mod_category", "nu", "sigma2", "alpha", habNames, CIlo, CIhi)
     
     # Create results dataframe
     nrow <- length(spp)
     ncol <- length(colNames)
     results <- data.frame(matrix(NA, nrow, ncol))
     colnames(results) <- colNames
     results$sp <- spp

     # Loop over each species, assign to category and estimate any parameters
     for(i in 1:length(spp)){
          
          # Write progress statement
          writeLines("###################")
          writeLines("###################")
          writeLines("###################")
          writeLines("###################")
          writeLines("")
          writeLines(paste("Run", i,"/", length(spp)))
          writeLines("")
          writeLines("###################")
          writeLines("###################")
          writeLines("###################")
          writeLines("###################")
          
          # Subset to a single species 
          sp <- dataframe[dataframe[,species_codes]==spp[i],]
          
          # Create point pattern object from stem map 
          stems <- as.ppp(cbind(sp[,x_coords], sp[,y_coords]), W=plotdim)
          
          ### CATEGORY 1 TEST: CSR ####
          
          # If there are no departures from CSR assign the species to C1
          # and move on to next iteration of the loop
          if(gof_test(stems, nsim=999)==TRUE){
               results[i,"mod_category"] <- "CSR"
               next 
          }
          
          ### CATEGORY 2 TEST: HABITAT ONLY ####
          
          # Specify the habitat trend formula 
          trend <- as.formula(paste("~", paste(names(hab_vars), collapse="+"), sep=""))
          
          # Fit heterogenious poisson model (HP)
          mod <- ppm(stems, trend=trend, covariates=hab_vars, control=list(maxit=10000))
          
          # Perform backwards selection via AIC, specifying AIC with k=2 
          step_mod <- step(mod, direction="backward", k=2) 
          
          # Check for deviations from the HP model
          hp_dev <- gof_test(step_mod, nsim=999)
          
          # If there are no significant departures from the HP model 
          # and the species is associated with at least one habitat attribute 
          if(hp_dev==TRUE & length(step_mod$coef)>1){
               
               # Assign it to the 'Habitat only' category
               results[i,"mod_category"] <- "H"
               
               # Add model coefficients to the results dataframe 
               results[i, names(step_mod$coef)[-1]] <- step_mod$coef[-1]     
               
               # Add confidence intervals to results dataframe
               modsum <- summary(step_mod)$coefs.SE.CI[-1,]
               results[which(results$sp==spp[i]), paste0(rownames(modsum), "CI95.lo")] <- modsum$CI95.lo
               results[which(results$sp==spp[i]), paste0(rownames(modsum), "CI95.hi")] <- modsum$CI95.hi
               
               # And move on to the next iteration of the loop
               next
               
          }
          
          ### CATEGORY 3 TEST: DISPERSAL ONLY ####
          # If there are significant departures from the HP model 
          # and the species is not associated with any habitat attribute 
          if(hp_dev==FALSE & length(step_mod$coef)==1){
               
               # Assign it to the 'Dispersal only' category
               results[i,"mod_category"] <- "D"
               
               # Fit a Thomas cluster process to the data
               mod <- kppm(stems ~ 1, "Thomas", statistic="pcf", method="mincon", 
                           model="matern", nu=nu, control=list(maxit=10000))
               
               # Add clustering paramters to results dataframe 
               results[i, c("nu", "sigma2", "alpha")] <- c(nu, mod$mu, mod$clustpar['scale'])
               
               # And move on to the next iteration of the loop
               next
               
          }
          
          ### CATEGORY 4 TEST: HABITAT AND DISPERSAL ####
          
          # If there are significant departures from the HP model 
          # and the species is associated at least one habitat attribute 
          if(hp_dev==FALSE & length(step_mod$coef) > 1){
               
               # Assign it to the 'Habitat and dispersal' category
               results[i,"mod_category"] <- "HD"
     
               # Extract names of significant parameters 
               sigPars <- names(step_mod$coef)[2:length(coef(step_mod))]
               
               # Specify the new habitat trend formula 
               trend <- as.formula(paste("~", paste(sigPars, collapse="+"), sep=""))
               
               # Remove insignificant habitat images 
               sig_hab_vars <- as.list(rep(NA, length(sigPars)))
               for(j in 1:length(sigPars)){
                    sig_hab_vars[[j]] <- hab_vars[[sigPars[j]]]
               }
               names(sig_hab_vars) <- sigPars
               
               # Fit Cox model with all habitat variables 
               mod <- kppm(stems, trend=trend, covariates=sig_hab_vars, "LGCP", statistic="pcf", method="mincon", 
                           model="matern", nu=nu, control=list(maxit=10000), improve.type="quasi")
               
               # Add model coefficients to the results dataframe 
               results[i, sigPars] <- as.numeric(coef(mod))[2:length(coef(mod))]   
               
               # Add confidence intervals to results dataframe
               modsum <- summary(mod)$coefs.SE.CI[-1,]
               results[i, paste0(rownames(modsum), "CI95.lo")] <- modsum$CI95.lo
               results[i, paste0(rownames(modsum), "CI95.hi")] <- modsum$CI95.hi
               
               # Add clustering paramters to results dataframe 
               results[i, c("nu", "sigma2", "alpha")] <- c(nu, mod$par['sigma2'], mod$par['alpha'])
               
          }
     }
     
     return(results)
     
}
