## Description: this simulation script is to estimate the proprotion of VUS that are likely pathogenic 
## (true_pi= 0.1 is specified at beginning). It also contains the code for calculating positive predictive value.


library(pracma)
library(maxLik)


script_num=1

# get script_num from command line
args = commandArgs(trailingOnly=T)
if (length(args)>=1) {
  script_num = as.numeric(args[1])
}

## Load simulated family pedigree data
load(paste0("../family_data/5k_families/5k_fam", script_num, ".RData"))

res.all = res.sub
non_carrier = res.all[which(res.all$BRCA1==0 & res.all$BRCA2==0), ]
carrier = res.all[which(res.all$BRCA1==1 | res.all$BRCA2==1),]

set.seed(script_num)

## Create reported genetic test results for all probands

# true_pi: the true proprotion of VUS that are likely pathogenic
# prop_VUS: proportion of individuals in the dataset with a reported genetic test result of VUS
# n_pos: number of positive probands to sample to be VUS probands 
# n_neg: number of negative probands to sample to be VUS probands 
# VUS: reported VUS probands 
# pos: reported positive probands 
# neg: reported negative probands 
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


## Logit transformation on carrier scores
x1 = log(pos$carrier_prob/(1-pos$carrier_prob)) 
x2 = log(neg$carrier_prob/(1-neg$carrier_prob))
x3 = log(VUS$carrier_prob/(1-VUS$carrier_prob)) 

## Get the range of carrier scores in reported VUS group
VUS_min = range(x3)[1]
VUS_max = range(x3)[2]


## Density function estimation 
df1 <- approxfun(density(x1, bw = "sj", from = VUS_min, to = VUS_max))
df2 <- approxfun(density(x2, bw = "sj", from = VUS_min, to = VUS_max))


## Likelihood approach to estimate 
lik <- function(pi) sum(log(pi*df1(x3)+ (1-pi)*df2(x3)))
A = matrix(c(1,-1),byrow = TRUE)
B = c(0,1)
ml <- maxLik(lik, start=0.4, constraints = list(ineqA=A, ineqB=B))
point_estimate = summary(ml)$estimate[1]

save(point_estimate, file=paste0("../results/5k_res_", script_num, ".RData"))


## Calculate positive predictive value
scores = seq(VUS_min, VUS_max, by = 0.01)
ppv = point_estimate*df1(scores)/(point_estimate*df1(scores) + (1-point_estimate)*df2(scores))

