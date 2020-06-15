
## Description: this script is to generate N family pedigree data using different family sizes
## and structures sampling from the USC-Standford HCP cohort. Carrier scores are calculated for all probands.

# load packages
library(BMmultigene)
library(BayesMendel)


script_num=1

# get script_num from command line
args = commandArgs(trailingOnly=T)
if (length(args)>=1) {
  script_num = as.numeric(args[1])
}

set.seed(script_num)


# load USC family structures. USC cohort contains 2000 families in total. 
load("usc_fam_structs.rData")
sample_structure_id = sample(1:2000, 1)

N = 1e5


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
  
# Bit of a hack: males don't get ovarian cancer, so there is no death from other causes
# column for that. Use death from other causes for female ovarian cancer. 
deathOtherCauses = death.othercauses[,c(2,1,3,3)] 
names(deathOtherCauses) = c("femaleBC", "maleBC", "femaleOC", "maleOC")
  
# Cancer
CP = genCancerPen(mutations, cancers, penCancersF, penCancersM)
  
# Other causes
ODP = genOtherDeathPen(cancers, deathOtherCauses)
  
  
### 
# Load data
data(compriskSurv)
  
# Separate competing risk by gender
comprisk = list(Female=compriskSurv[,1:4], Male=compriskSurv[,5:8])
  
# Rename columns to match possible genotypes
colnames(comprisk$Female) = CP$PG 
colnames(comprisk$Male) = CP$PG
  
  
# Specify number of males and females in each branch of the family. 
nSibsPatern = usc_fam_structs[[sample_structure_id]][1, ]  # paternal female/male relatives
nSibsMatern = usc_fam_structs[[sample_structure_id]][2, ]  # maternal female/male relatives
nSibs = usc_fam_structs[[sample_structure_id]][3, ] # number of siblings
nGrandchild = usc_fam_structs[[sample_structure_id]][4, ]
  
  
## simulate families
fam = sim.simFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, prevs, CP, includeGeno=TRUE)
lik = estLik(fam, CP, ODP, comprisk)
res= data.frame(pp.peelingParing(fam, prevs, lik, 2)) 
  
  
for (i in 2:N){
      
   sample_structure_id = sample(1:2000, 1)
   
   # Specify number of males and females in each branch of the family.
   nSibsPatern = usc_fam_structs[[sample_structure_id]][1, ]  # paternal female/male relatives
   nSibsMatern = usc_fam_structs[[sample_structure_id]][2, ]  # maternal female/male relatives
   nSibs = usc_fam_structs[[sample_structure_id]][3, ] # number of siblings
   nGrandchild = usc_fam_structs[[sample_structure_id]][4, ]
   
   fam_new = sim.simFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, prevs, CP, includeGeno=TRUE)
   lik = estLik(fam_new, CP, ODP, comprisk)
   res_new= data.frame(pp.peelingParing(fam_new, prevs, lik, 2))
    
   fam = rbind(fam,fam_new)
   res = cbind(res,res_new)
  }
  

# re-format simulated family data and store it. 
res2 = t(res)
carrier_prob = 1-res2[,1]
proband = fam[which(fam$isProband==1),]
fam.id = c(1:N)
res.all = cbind(fam.id, proband, carrier_prob)
res.all = res.all[, c('fam.id','BRCA1','BRCA2','carrier_prob')]

save(res.sub, file=paste0("../family_data/results/1k_fam", script_num, ".RData"))

  
