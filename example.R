#####################
#
# Example of a single run of the model
# Dominique Gravel
# June 2nd, 2016
# 
####################
rm(list = ls())

library(metadym)
source("parameters.R")

# Run the model
nsteps = 1000
set.seed(1)

# running without SA

out = main(XY,E,Y0,pars,A,nsteps)

computeAlphaGamma(out)
plotOutput(out)

pars$SA = 0.5

out = main(XY,E,Y0,pars,A,nsteps)

computeAlphaGamma(out)
plotOutput(out)





plot(XY[,1], XY[,2], cex = 2*E[,1])
library(ape)


ozone.dists <- as.matrix(dist(XY))

ozone.dists.inv <- 1/ozone.dists
diag(ozone.dists.inv) <- 0
Moran.I(E[,1], ozone.dists.inv)
