
## Script to generate Figure 1.
getwd()

library(BMmultigene)
library(BayesMendel)
library(pracma)
library(maxLik)

# load a simulated cohort of sample size = 5000
load("../5k_fam1.RData")

res.all = res.sub
non_carrier = res.all[which(res.all$BRCA1==0 & res.all$BRCA2==0), ]
carrier = res.all[which(res.all$BRCA1==1 | res.all$BRCA2==1),]


## Code for plotting Figure 1. 
par(mfrow = c(2,2))

## Carrier score distribution when theta = 0 
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

# Density function estimation on original scale
x1 = pos$carrier_prob 
x2 = neg$carrier_prob 
x3 = VUS$carrier_prob  

df1 <- approxfun(density(x1, bw = "sj", from = 0, to =1))
df2 <- approxfun(density(x2, bw = "sj", from = 0, to =1))

plot(density(x2, bw = "sj"), col = 2,
     main = "", xlab = expression(paste("(a) ", theta, '=0')) ) #black
lines(density(x1, bw = "sj"), col = 1) #red
lines(density(x3, bw = "sj"), col =3)  #green
legend("topright",legend = c("Positive", "Negative","VUS"), col = c(1:3), pch = 20)

# Density function estimation on logit scale
lx1 = log(pos$carrier_prob/(1-pos$carrier_prob)) 
lx2 = log(neg$carrier_prob/(1-neg$carrier_prob)) 
lx3 = log(VUS$carrier_prob/(1-VUS$carrier_prob)) 

logit_df1 <- approxfun(density(lx1, bw = "sj", from = -8, to =11))
logit_df2 <- approxfun(density(lx2, bw = "sj", from = -8, to =11))

plot(density(lx1, bw = "sj"), col = 1, ylim = c(0,0.6),
     main = "", xlab = expression(paste("(b) ", theta, '=0')) ) #black
lines(density(lx2, bw = "sj"), col = 2) #red
lines(density(lx3, bw = "sj"), col =3)  #green
legend("topright",legend = c("Positive", "Negative","VUS"), col = c(1:3), pch = 20)



## Carrier score distribution when theta = 0.05
true_pi = 0.1
prop_VUS = 0.3
theta = 0.05

n_pos = nrow(res.all)*prop_VUS*true_pi
n_neg = nrow(res.all)*prop_VUS*(1-true_pi)

famid1 = sample (carrier$fam.id, n_pos, replace = FALSE)
famid2 = sample (non_carrier$fam.id, n_neg, replace = FALSE)

index1 = which(carrier$fam.id %in% famid1)
index2 =  which(non_carrier$fam.id %in% famid2)

VUS = rbind(carrier[index1,], non_carrier[index2,])
pos = carrier[-index1,]
neg = non_carrier[-index2,]


n_falneg = round(nrow(neg)*theta/(1-theta))
famid3 = sample(pos$fam.id, n_falneg, replace = FALSE)
index3 = which(pos$fam.id %in% famid3)

pos2 = pos[-index3,]
neg2 = rbind(neg, pos[index3,])


# Density function estimation on original scale
x1 = pos$carrier_prob 
x2 = neg$carrier_prob 
x3 = VUS$carrier_prob  

df1 <- approxfun(density(x1, bw = "sj", from = 0, to =1))
df2 <- approxfun(density(x2, bw = "sj", from = 0, to =1))

plot(density(x2, bw = "sj"), col = 2,
     main = "", xlab = expression(paste("(c) ", theta, '=0.05'))   ) #black
lines(density(x1, bw = "sj"), col = 1) #red
lines(density(x3, bw = "sj"), col =3)  #green
legend("topright",legend = c("Positive", "Negative","VUS"), col = c(1:3), pch = 20)


# Density function estimation on logit scale
lx1 = log(pos$carrier_prob/(1-pos$carrier_prob)) 
lx2 = log(neg$carrier_prob/(1-neg$carrier_prob)) 
lx3 = log(VUS$carrier_prob/(1-VUS$carrier_prob)) 

logit_df1 <- approxfun(density(lx1, bw = "sj", from = -8, to =11))
logit_df2 <- approxfun(density(lx2, bw = "sj", from = -8, to =11))

plot(density(lx1, bw = "sj"), col = 1, ylim = c(0,0.6),
     main = "", xlab = expression(paste("(d) ", theta, '=0.05')) ) #black
lines(density(lx2, bw = "sj"), col = 2) #red
lines(density(lx3, bw = "sj"), col =3)  #green
legend("topright",legend = c("Positive", "Negative","VUS"), col = c(1:3), pch = 20)
