
## Ian McFadden
## imcfadden@ucla.edu

################################################################################
#
# gof_test
#
# Tests for departures from complete spatial randomness (CSR) or departures from
# a fitted model using the Loosmore's goodness of fit test 
#
# If observed pattern is indistinguishable from CSR or the supplied model
# the p-value is above 0.05 and function returns TRUE indicating there are
# no departures  
#
# If the p-value is 0.05 or below the function returns FALSE indicating that 
# there are significant departures from either CSR or the supplied model 
#
# The function uses three different inter-point correlation metrics: 
#    pcf- the pair correlation function 
#    Lest- the L function 
#    Fest- the empty space function 
#
# ARGUMENTS
# obj     Either the spatial point pattern in the case of tests for departures 
#         from CSR or a fitted spatial point pattern model
# nsim    Number of simulations to run for the null distribution 
#
################################################################################

gof_test <- function(obj, nsim=999){
     
     require(spatstat)
     
     # Make results vector
     res <- rep(NA,3)
     names(res) <- c("pcf", "L_func", "empty_space")
     
     # Perfom tests and fill results vector 
     res["pcf"] <- dclf.test(obj, fun="pcf", nsim=nsim)$p.value 
     res["L_func"] <- dclf.test(obj, fun="Lest", nsim=nsim)$p.value  
     res["empty_space"] <- dclf.test(obj, fun="Fest", nsim=nsim)$p.value 
     
     # Correct for mutliple comparisons 
     res <- p.adjust(res, "BH")
     
     # Return the output 
     if(min(res) > 0.05){
          TRUE
     } else {
          FALSE
     }
}
