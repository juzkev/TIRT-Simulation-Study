##----------------------------------------------------------------------------------
##
##  Author:       Kevin Koh
##  Description:  Runs simulation on high correlation
##
##----------------------------------------------------------------------------------

library(MplusAutomation)

# Create Monte Carlo template----
tmp <-
  '[[init]]
iterators = model;
model = 1:7;
corr#model = 0.3 0.4 0.5 0.6 0.7 0.8 0.9;
filename = "Monte Carlo 4 traits SIMULATION high correlation [[corr#model]].inp";
outputDirectory = Data/highcorr/[[corr#model]];
[[/init]]

TITLE: Generates 4 traits with high correlation of [[corr#model]]
Generates enough errors for 12 items with error=1 (e) and 12 items with error=0.5 (se)
MONTECARLO:

NAMES = trait1-trait4 e1-e12 se1-se12;

NOBSERVATIONS = 2000;
NREPS = 1000;
SEED = 2204;
REPSAVE=ALL;
SAVE = 4TraitsHighCorr[[corr#model]]Rep*.dat;

MODEL POPULATION:

!standardised factors
trait1-trait4@1;

!factors are uncorrelated with errors, and errors are uncorrelated with each other
!factors are correlated
trait1 WITH trait2@-0.4 trait3@0 trait4@.4;
trait2 WITH trait3@[[corr#model]] trait4@-.3;
trait3 WITH trait4@0;

!Errors
e1-e12@1;
se1-se12@0.5;'

## Write to file----
write(tmp, "highcorrtmplt.txt")
createModels("highcorrtmplt.txt")

# Create analysis template----
tmp <-
  '[[init]]
iterators = model;
model = 1:7;
corr#model = 0.3 0.4 0.5 0.6 0.7 0.8 0.9;
filename = "3Traits4Triplets high correlation [[corr#model]].inp";
outputDirectory = Data/highcorr/[[corr#model]];
[[/init]]


TITLE: Model with 4 triplets measuring 3 traits at high correlation of [[corr#model]]
DATA:   !FILE IS 4TraitsRep5.dat;
FILE IS 4TraitsHighCorr[[corr#model]]Replist.dat; TYPE=MONTECARLO;

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
write(tmp, "highcorranatmplt.txt")
createModels("highcorranatmplt.txt")

# Run data generation----
runModels(getwd(), recursive = TRUE, filefilter = "Monte Carlo 4 traits SIMULATION high correlation \\d\\.\\d\\.inp")

# Run analysis----
runModels(getwd(), recursive = TRUE, filefilter = "3Traits4Triplets high correlation \\d\\.\\d\\.inp")

# Read output----
allOutput <-
  readModels(getwd(), recursive = TRUE, filefilter = "3traits4triplets high correlation \\d\\.\\d")

# Compile output into a table in parameters----
# Select only the needed fit indexes
modelSum = data.frame(matrix(ncol = 16, nrow = 7))
names(modelSum) <- c(
  'Parameters',
  'ChiSqM_DF',
  'ChiSqM_Mean',
  'ChiSqM_SD',
  'ChiSqM_NumComputations',
  'RMSEA_Mean',
  'RMSEA_SD',
  'RMSEA_NumComputations',
  'meanLoadingsRelativeBias',
  'seLoadingsRelativeBias',
  'meanCorrelationsAbsoluteBias',
  'seCorreleationsAbsoluteBias',
  'meanUniquenessRelativeBias',
  'seUniquenessRelativeBias',
  'meanThresholdsRelativeBias',
  'seThresholdsRelativeBias'
)
row.names(modelSum) <- paste0("Model", 1:7)
modelParam <- list()
tmpregex <- c(".*BY", ".*WITH", ".+Variances", "Thresholds")
popThrshlds <-
  c(0.5, -1.2, -1.7, 0.7, 1, 0.3,-0.7,-1.2,-0.5, 0.7, 1.2, 0.5)

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
  # calculate relative bias for loadings and uniqueness
  modelParam[[i]]$parambias <-
    abs(modelParam[[i]]$average - modelParam[[i]]$population) / abs(modelParam[[i]]$population)
  modelParam[[i]]$sebias <-
    abs(modelParam[[i]]$average_se - modelParam[[i]]$population_sd) / abs(modelParam[[i]]$population_sd)
  for (j in 1:4) {
    modelSum[i, names(modelSum)[j*2 - 1 + 8]] <-
      mean(modelParam[[i]][with(modelParam[[i]], grepl(tmpregex[j], paramHeader)) , 'parambias'])
    modelSum[i, names(modelSum)[j*2 + 8]] <-
      mean(modelParam[[i]][with(modelParam[[i]], grepl(tmpregex[j], paramHeader)) , 'sebias'])
  }
  # Population value for Thresholds not correctly specified in Mplus output, thus using manually defined threshold
  modelSum[i, 'meanThresholdsRelativeBias'] <-
    mean(abs(modelParam[[i]][with(modelParam[[i]], grepl("Thresholds", paramHeader)) , 'average'] - popThrshlds) / abs(popThrshlds))
  # calculate absolute bias for correlations, as some correlations are 0
  modelSum[i, 'meanCorrelationsAbsoluteBias'] <-
    mean(abs(modelParam[[i]][with(modelParam[[i]], grepl(".*WITH", paramHeader)) , 'average'] - modelParam[[i]][with(modelParam[[i]], grepl(".*WITH", paramHeader)) , 'population']))
  # calculate absolute SE bias for correlations, as some correlations are 0
  modelSum[i, 'seCorreleationsAbsoluteBias'] <-
    mean(abs(modelParam[[i]][with(modelParam[[i]], grepl(".*WITH", paramHeader)) , 'average_se'] - modelParam[[i]][with(modelParam[[i]], grepl(".*WITH", paramHeader)) , 'population_sd']))
  # adjust for triplets redundencies r=n(n-1)(n-2)/6 = 1 per block
  modelSum[i, 'ChiSqM_DFadj'] <-
    (modelSum[i, 'ChiSqM_DF'] - 4) 
  # calculate Chi-sq probablity based on new DF
  modelSum[i, 'ChisqM_p'] <-
    (pchisq(modelSum[i,'ChiSqM_Mean'],modelSum[i, 'ChiSqM_DFadj'],lower.tail = FALSE)*2)
}

# Plot visualisation, Idk what to visualise sia----
(modelSum)
write.csv(modelSum, file = "High Corr.csv")
