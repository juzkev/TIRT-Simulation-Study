##----------------------------------------------------------------------------------
##  
##  Author:       Kevin Koh
##  Description:  Runs simulation on non-normal factor
##
##----------------------------------------------------------------------------------

cond <- c(0, 0.4, 0.8, 1.6, 2) # skewness kurtosis values from (Blanca, Arnau, López-Montiel, Bono, Bendayan, 2013)

library('simsem')

path.BE <- matrix(0, 28, 28)
BE <- bind(path.BE)

cov.RPS <- matrix(0, 28, 28)
cov.RPS[1, 2:4] <- cov.RPS[2:4, 1] <- c(NA, 0, NA)
cov.RPS[2, 3:4] <- cov.RPS[3:4, 2] <-  NA
val.RPS <- diag(1, 4)
val.RPS[1, 2] <- val.RPS[2, 1] <- -0.4
val.RPS[1, 4] <- val.RPS[4, 1] <- 0.4
val.RPS[2, 3] <- val.RPS[3, 2] <- 0.3
val.RPS[2, 4] <- val.RPS[4, 2] <- -0.3
valbig.RPS <- diag(1, 28)
valbig.RPS[1:4, 1:4] <- val.RPS
valbig.RPS[17:28, 17:28] <- diag(0.5, 12)
RPS <- binds(cov.RPS, valbig.RPS)

path.Model <- model(BE = BE,
                    RPS = RPS,
                    modelType = "Path")

n1 <- list(mean = 0, sd = 1)
n2 <- list(mean = 0, sd = 0.5)

wd <- getwd()
dir.create(file.path(wd, 'Data', 'nnorm'),
           showWarnings = FALSE,
           recursive = TRUE)

for (i in c(1:length(cond))) {
  cat('Monte Carlo data generated has started for model ', i, '.\n', sep = "")
  skewness <- c(rep(cond[i], 3), rep(0, (25))) # insert i here
  kurtosis <- c(rep(cond[i], 3), rep(0, (25))) # insert i here
  facDist <-
    bindDist(skewness = skewness,
             kurtosis = kurtosis)
  tmp.dir <- file.path(wd, 'Data', 'nnorm', paste0("model ", i))
  dir.create(tmp.dir, showWarnings = FALSE, recursive = TRUE)
  setwd(tmp.dir)
  exportData(
    nRep = 1000,
    model = path.Model,
    n = 2000,
    program = "Mplus",
    fileStem = paste0("4Traitsnnorm", i, "Rep") ,
    sequential = TRUE,
    facDist = facDist,
    seed = 2204
  )
  cat('Monte Carlo data generated completed for model ', i, '.\n', sep = "")
}
setwd(wd)

library(MplusAutomation)

# Create analysis template----
tmp <- 
  '[[init]]
iterators = model;
model = 1:5;
filename = "3Traits4Triplets nnorm [[model]].inp";
outputDirectory = "Data/nnorm/model [[model]]";
[[/init]]


TITLE: Model with 4 triplets measuring 3 traits with non-normal distrubution, skewness kurtosis of [[model]]
DATA:   !FILE IS 4TraitsRep5.dat;
FILE IS 4Traitsnnorm[[model]]Rep.dat; TYPE=MONTECARLO;

VARIABLE:
NAMES = trait1-trait4 e1-e12 se1-se12;
USEVARIABLES=
i1i2 i1i3 i2i3
i4i5 i4i6 i5i6
i7i8 i7i9 i8i9
i10i11 i10i12 i11i12;
CATEGORICAL = i1i2-i11i12;

DEFINE:

i1i2=   -0.5+   (1)*trait1-     (0.8)*trait2+e1-e2;
i1i3=   1.2+    (1)*trait1-     (1.3)*trait3+e1-e3;
i2i3=   1.7+    (0.8)*trait2-   (1.3)*trait3+e2-e3;
i4i5=   -0.7+   (-1.3)*trait1-  (1)*trait2+e4-e5;
i4i6=   -1+     (-1.3)*trait1-  (0.8)*trait3+e4-e6;
i5i6=   -0.3+   (1)*trait2-     (0.8)*trait3+e5-e6;
i7i8=   0.7+    (0.8)*trait1-   (1.3)*trait2+e7-e8;
i7i9=   1.2+    (0.8)*trait1-   (-1)*trait3+e7-e9;
i8i9=   0.5+    (1.3)*trait2-   (-1)*trait3+e8-e9;
i10i11= -0.7+   (1.3)*trait1-   (-0.8)*trait2+e10-e11;
i10i12= -1.2+   (1.3)*trait1-   (1)*trait3+e10-e12;
i11i12= -0.5+   (-0.8)*trait2-  (1)*trait3+e11-e12;


CUT i1i2-i11i12 (0);

ANALYSIS:
ESTIMATOR=ulsmv;
PARAMETERIZATION=THETA;

MODEL:

Trait1  BY  	
i1i2*1  (L1)
i1i3*1  (L1)
i4i5*-1.3  (L4)
i4i6*-1.3  (L4)
i7i8*.8  (L7)
i7i9*.8  (L7)
i10i11*1.3  (L10)
i10i12*1.3  (L10);
Trait2  BY  	
i1i2*-.8  (L2_n)
i2i3*.8  (L2)
i4i5*-1  (L5_n)
i5i6*1  (L5)
i7i8*-1.3  (L8_n)
i8i9*1.3  (L8)
i10i11*.8  (L11_n)
i11i12*-.8  (L11);
Trait3  BY  	
i1i3*-1.3  (L3_n)
i2i3*-1.3  (L3_n)
i4i6*-.8  (L6_n)
i5i6*-.8  (L6_n)
i7i9*1  (L9_n)
i8i9*1  (L9_n)
i10i12*-1  (L12_n)
i11i12*-1  (L12_n);

! variances for all traits are set to 1	
Trait1-Trait3@1;	

! starting values for correlations between traits	
Trait1 WITH Trait2*-0.8	Trait3*0;
Trait2 WITH Trait3*0.3;

! declare uniquenesses 	
i1i2*2 (e1e2);	
i1i3*2 (e1e3);	
i2i3*2 (e2e3);	
i4i5*2 (e4e5);	
i4i6*2 (e4e6);	
i5i6*2 (e5e6);	
i7i8*2 (e7e8);	
i7i9*2 (e7e9);	
i8i9*2 (e8e9);
i10i11*2 (e10e11);
i10i12*2 (e10e12);
i11i12*2 (e11e12);	

! declare correlated uniqunesses and set their starting values	
i1i2 WITH i1i3*1 (e1);	
i1i2 WITH i2i3*-1 (e2_n);	
i1i3 WITH i2i3*1 (e3);	

i4i5 WITH i4i6*1 (e4);	
i4i5 WITH i5i6*-1 (e5_n);	
i4i6 WITH i5i6*1 (e6);	

i7i8 WITH i7i9*1 (e7);	
i7i8 WITH i8i9*-1 (e8_n);	
i7i9 WITH i8i9*1 (e9);	

i10i11 WITH i10i12*1 (e10);	
i10i11 WITH i11i12*-1 (e11_n);	
i10i12 WITH i11i12*1 (e12);	

MODEL CONSTRAINT:	

!factor loadings relating to the same item are equal in absolute value	
L2_n = -L2;	
L5_n = -L5;	
L8_n = -L8;	
L11_n = -L11;

! pair uniqueness is equal to sum of 2 utility uniqunesses	
e1e2 = e1 - e2_n;	
e1e3 = e1 + e3;	
e2e3 = -e2_n + e3;	
e4e5 = e4 - e5_n;	
e4e6 = e4 + e6;	
e5e6 = -e5_n + e6;	
e7e8 = e7 - e8_n;	
e7e9 = e7 + e9;	
e8e9 = -e8_n + e9;	
e10e11 = e10 - e11_n;
e10e12 = e10 + e12;
e11e12 = -e11_n + e12;	

! fix one uniqueness per block for identification	
e1=1;
e4=1;
e7=1;
e10=1;

OUTPUT: STDY;

!SAVEDATA: file is triplets.dat;
!SAVE=FSCORES;'

# Write to file----
write(tmp, "nnormanatmplt.txt")
createModels("nnormanatmplt.txt")

# Run analysis----
runModels(getwd(), recursive = TRUE, filefilter = "3Traits4Triplets nnorm \\d\\.inp")

# Read output----
allOutput <-
  readModels(getwd(), recursive = TRUE, filefilter = "3Traits4Triplets nnorm \\d")

# Compile output into a table in parameters----
# Select only the needed fit indexes
modelSum = data.frame(matrix(ncol = 12, nrow = length(cond)))
names(modelSum) <- c(
  'Parameters',
  'ChiSqM_DF',
  'ChiSqM_Mean',
  'ChiSqM_SD',
  'ChiSqM_NumComputations',
  'RMSEA_Mean',
  'RMSEA_SD',
  'RMSEA_NumComputations',
  'meanLoadingsBias',
  'meanCorrelationsBias',
  'meanUniquenessBias',
  'meanThresholdsBias'
)
row.names(modelSum) <- paste0("Model", 1:length(cond))
modelParam <- list()
tmpregex <- c(".*BY", ".*WITH", ".+Variances", "Thresholds")
popThrshlds <-
  c(0.5,-1.2,-1.7, 0.7, 1, 0.3, -0.7, -1.2, -0.5, 0.7, 1.2, 0.5)

# Loop through each model
for (i in 1:(length(allOutput))) {
  modelSum[i, names(modelSum)[1:8]] <-
    allOutput[[i]]$summaries[names(modelSum)[1:8]]
  modelParam[[i]] <-
    allOutput[[i]]$parameters$unstandardized[c('paramHeader',
                                               'param',
                                               'population',
                                               'average',
                                               'population_sd',
                                               'average_se')]
  # calculate bias
  modelParam[[i]]$bias <-
    abs(modelParam[[i]]$average - modelParam[[i]]$population)
  for (j in seq_along(names(modelSum)[9:12])) {
    modelSum[i, names(modelSum)[j + 8]] <-
      mean(modelParam[[i]][with(modelParam[[i]], grepl(tmpregex[j], paramHeader)) , 'bias'])
  }
  modelSum[i, 'meanThresholdsBias'] <-
    mean(abs(modelParam[[i]][with(modelParam[[i]], grepl("Thresholds", paramHeader)) , 'average'] - popThrshlds))
  modelSum[i, 'ChiSqM_DFadj'] <-
    (modelSum[i, 'ChiSqM_DF'] - 4) # adjust for triplets redundencies r=n(n-1)(n-2)/6 = 1 per block
}

# Plot visualisation, Idk what to visualise sia----
(modelSum)
write.csv(modelSum, file = "nnorm.csv")
