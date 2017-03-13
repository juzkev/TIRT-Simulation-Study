##----------------------------------------------------------------------------------
##  
##  Author:       Kevin Koh
##  Description:  Generating normal data for simulation studies according to the
##                parameters specified for Brown(2012). Variables are generated in 
##                the following order and specification: n1-4, e1-12, se1-12.
##                n1~~n2=-0.4
##                  ~~n4=0.4
##                n2~~n3=0.3
##                  ~~n4=0.3
##  Version:      1.0
##
##----------------------------------------------------------------------------------

library('semPlot')

## Load the library for SEM simulation
library('simsem')

# Regression coefficients (1 directional arrow) between endogenous factors (IVs)
# all regression coefficient 
path.BE <- matrix(0, 28, 28)
BE <- bind(path.BE)

## Fixing (with value) and freeing (NA) parameters, freed correlations between the n factors
cov.RPS <- matrix(0, 28, 28)
cov.RPS[1, 2:4] <- cov.RPS[2:4, 1] <- c(-0.4, 0, 0.4)  #c(NA, 0, NA) ## not sure to free it or to specify
cov.RPS[2, 3:4] <- cov.RPS[3:4, 2] <- c(0.3, 0.3) #NA ## not sure to free it or to specify

## Residual correlation matrix amoung endogenous factors (n1 to n4)
val.RPS <- diag(1, 4) # variance 1
val.RPS[1, 2] <- val.RPS[2, 1] <- -0.4 # n1~~n2=-0.4
val.RPS[1, 4] <- val.RPS[4, 1] <- 0.4 # n1~~n4=0.4
val.RPS[2, 3] <- val.RPS[3, 2] <- 0.3 # n2~~n3=0.3
val.RPS[2, 4] <- val.RPS[4, 2] <- -0.3 # n2~~n4=0.3
val.RPS # check matrix

## Residual correlation matrix amoung endogenous factors (n1 to n4, e1-12, se1-12)
valbig.RPS <- matrix(0, 28, 28)
valbig.RPS[1:4, 1:4] <- val.RPS
valbig.RPS[5:16, 5:16] <- diag(1, 12) # e1-12 variance 1
valbig.RPS[17:28, 17:28] <- diag(0.5, 12) # se1-12 variance 0.5
valbig.RPS # check matrix
RPS <- binds(free = cov.RPS, popParam = valbig.RPS)

path.Model <- model(BE = BE, RPS = RPS, modelType = "Path")

semPaths(path.Model@pt, what = 'est', intercepts = FALSE, layout = 'tree', residuals = FALSE, curvature = 2 , curve = 2) ## values does not seem to come out right if i don't specify it in cov.RPS

# specify the distribution for the e1-12 and se1-12
n1 <- list(mean = 0, sd = 1)
n2 <- list(mean = 0, sd = 0.5)

# Input distribution type
facDist <-
  bindDist(
    c(rep("norm", 28)),
    n1,
    n1,
    n1,
    n1,
    n1,
    n1,
    n1,
    n1,
    n1,
    n1,
    n1,
    n1,
    n1,
    n1,
    n1,
    n1,
    n2,
    n2,
    n2,
    n2,
    n2,
    n2,
    n2,
    n2,
    n2,
    n2,
    n2,
    n2
  )

exportData(
  nRep = 100,
  model = path.Model,
  n = 2000,
  program = "Mplus",
  fileStem = "Norm4TraitsRep",
  sequential = TRUE,
  facDist = facDist,
  seed = 2204
)

