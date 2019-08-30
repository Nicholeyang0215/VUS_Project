

getwd()

library(BMmultigene)
library(BayesMendel)
library(pracma)
library(maxLik)


seed <- round(runif(1, 0, 5000))
set.seed(seed)



  ####################################
  ####### Jane's simulation code #####
  ####################################
  
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
  
  #####################################
  ######## Simulate 5000 families #####
  #####################################
  
  fam = sim.simFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, prevs, CP, includeGeno=TRUE)
  lik = estLik(fam, CP, ODP, comprisk)
  res= data.frame(pp.peelingParing(fam, prevs, lik, 2)) 
  
  
  for (i in 2:1000){
    
    fam_new = sim.simFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, prevs, CP, includeGeno=TRUE)
    lik = estLik(fam_new, CP, ODP, comprisk)
    res_new= data.frame(pp.peelingParing(fam_new, prevs, lik, 2))
    
    fam = rbind(fam,fam_new)
    res = cbind(res,res_new)
  }
  
  res2 = t(res)
  carrier_prob = 1-res2[,1]
  proband = fam[which(fam$isProband==1),]
  
  
  ############
  fam.id = c(1:1000)
  res.all = cbind(fam.id,proband,carrier_prob)
  non_carrier = res.all[which(res.all$BRCA1==0 & res.all$BRCA2==0), ]
  carrier = res.all[which(res.all$BRCA1==1 | res.all$BRCA2==1),]

  
  ################ Create VUS, positive, and negative #########
  true_pi = 0.1
  prop_VUS = 0.3
  
  n_pos = nrow(res.all)*prop_VUS*true_pi
  n_neg = nrow(res.all)*prop_VUS*(1-true_pi)
  
  famid1 = sample (carrier$fam.id, n_pos, replace = FALSE)
  famid2 = sample (non_carrier$fam.id, n_neg, replace = FALSE)
  
  index1 = which(carrier$fam.id %in% famid1)
  index2 =  which(non_carrier$fam.id %in% famid2)
  
  VUS = rbind(carrier[index1,], non_carrier[index2,])
  pos = carrier[-index1,]
  neg = non_carrier[-index2,]
  
  
  #### observed ones: VUS, pos2, neg2 
  theta = 0.05
  n_falneg = round(nrow(neg)*theta/(1-theta))
  famid3 = sample(pos$fam.id, n_falneg, replace = FALSE)
  index3 = which(pos$fam.id %in% famid3)
  
  pos2 = pos[-index3,]
  neg2 = rbind(neg, pos[index3,])

  
  ############### Density function estimation ##########
  x1 = log(pos2$carrier_prob/(1-pos2$carrier_prob)) ###positive carrier probabilities 
  x2 = log(neg2$carrier_prob/(1-neg2$carrier_prob)) ###negative carrier probabilities 
  x3 = log(VUS$carrier_prob/(1-VUS$carrier_prob)) 
  
  df1 <- approxfun(density(x1, bw = "sj", from = -8, to =11))
  df2 <- approxfun(density(x2, bw = "sj", from = -8, to =11))
  
  
  ########## maximum likelihood approach #################
  theta =0.05
  
  lik <- function(pi) sum(log(pi*df1(x3)+(1-pi)*((df2(x3) - theta*df1(x3))/(1-theta))))
  
  A = matrix(c(1,-1),byrow = TRUE)
  B = c(0,1)
  ml <- maxLik(lik, start=0.7, constraints = list(ineqA=A, ineqB=B))
  point_estimate1 = summary(ml)$estimate[1]

save(point_estimate1, file = paste("results_1k",seed,".RData", sep=""))

