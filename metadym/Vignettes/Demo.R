## ---- echo = F-----------------------------------------------------------
set.seed(1)

## ---- eval = F-----------------------------------------------------------
#  library(devtools)
#  install_url("https://dl.dropboxusercontent.com/s/cvsvuial6nqmqov/metadym_0.0.0.9000.tar.gz", dependencies = T)

## ------------------------------------------------------------------------
library(metadym)
?metadym

## ------------------------------------------------------------------------

# Global model parameters

nsteps = 500 # Number of time steps 
N = 200 # Number of patches
D = 2 # Number of niches
R = 25 # Number of species

# Process parameters

pars = getDefaultParameters() # type getDefaultParameters in console to see how defaults are created

# Changing process defaults

pars$SA = 0.04


