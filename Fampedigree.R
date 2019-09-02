

# Script to generate 5000 family pedigree data and calculate carrier scores for all probands. 

getwd()
library(BMmultigene)
library(BayesMendel)


seed <- round(runif(1, 0, 5000))
set.seed(seed)


# Consider breast cancer and ovarian cancer
cancers = c("BC", "OC")

# Consider BRCA1 and BRCA2 mutations
mutations = c("BRCA1", "BRCA2")

# BRCA1 and BRCA2 prevalences 
prevs = c(BRCA1=0.05, BRCA2=0.05) 

# Load data
data(BRCApenet.metaDSL.2008) 

# Get penetrance corresponding to noncarriers, BRCA1 carriers, and BRCA2 carriers
# Females 
penCancersF = list(BC=BRCApenet.metaDSL.2008$fFX[,c(1,2,4)], # Breast cancer
                   OC=BRCApenet.metaDSL.2008$fFY[,c(1,2,4)]) # Ovarian cancer
# Males
penCancersM = list(BC=BRCApenet.metaDSL.2008$fMX[,c(1,2,4)], # Breast cancer
                   OC=BRCApenet.metaDSL.2008$fMY[,c(1,2,4)]) # Ovarian cancer

data(death.othercauses) 

# Males don't get ovarian cancer, so there is no death from other causes
# Use death from other causes for female ovarian cancer. 
deathOtherCauses = death.othercauses[,c(2,1,3,3)] 
names(deathOtherCauses) = c("femaleBC", "maleBC", "femaleOC", "maleOC")

# Cancer
CP = genCancerPen(mutations, cancers, penCancersF, penCancersM)

# Other causes
ODP = genOtherDeathPen(cancers, deathOtherCauses)

# Load data
data(compriskSurv)

# Separate competing risk by gender
comprisk = list(Female=compriskSurv[,1:4], Male=compriskSurv[,5:8])

# Rename columns to match possible genotypes
colnames(comprisk$Female) = CP$PG 
colnames(comprisk$Male) = CP$PG


# Specify number of males and females in each branch of the family. 
nSibsPatern = c(1,2) # One paternal aunt and two paternal uncles 
nSibsMatern = c(2,2) # Two maternal aunts and two maternal uncles
nSibs = c(2,3) # Proband has two sisters and three brothers 
nGrandchild = c(2,1) # Each person in proband's generation has two daughters and a son 


# Generation of one family pedigree with carrier scores.
# Function sim.simFam generates a completed family matrix 
fam = sim.simFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, prevs, CP, includeGeno=TRUE)
# Function eskLik estimates the likelihood matrix, in other words, calculate probands carrier probability. 
lik = estLik(fam, CP, ODP, comprisk)
res= data.frame(pp.peelingParing(fam, prevs, lik, 2)) 


# The number of family pedigree to generate 
n = 5000

# Loop to combine the rest of family pedigrees 
for (i in 2:n){
  
  fam_new = sim.simFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, prevs, CP, includeGeno=TRUE)
  lik = estLik(fam_new, CP, ODP, comprisk)
  res_new= data.frame(pp.peelingParing(fam_new, prevs, lik, 2))
  
  fam = rbind(fam,fam_new)
  res = cbind(res,res_new)
}

res2 = t(res)
carrier_prob = 1-res2[,1]
proband = fam[which(fam$isProband==1),]


fam.id = c(1:n)
res.all = cbind(fam.id,proband,carrier_prob)
non_carrier = res.all[which(res.all$BRCA1==0 & res.all$BRCA2==0), ]
carrier = res.all[which(res.all$BRCA1==1 | res.all$BRCA2==1),]

