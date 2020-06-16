
## Script to visualize simulation results of positive predictive value (PPV). 


# load PPVs for all simulated cohorts
setwd("/Users/nicholeyang/Desktop/research2/VUSpaper_Revision/new_sim_results/Fig4/results")
res.list = list.files(pattern="ppv*")

n = 2000 ## number of replicates

# create carrier scores for the same length. 
sequence = seq(-9, 10, by = 0.01)

# create matrix to store postivie predictive value(dt_ppv). Columns represent replicates.
dt_ppv = matrix(NA, ncol = n, nrow = length(sequence)) 

## fill ppv into dt_ppv
for (i in 1:n){
  
  load(res.list[i])

  scores2 = round(scores,2)  ## round simulated scores to 2-digit
  indx = which(sequence == scores2[1])   ## find the place to start 
  
  for (j in c(1:length(scores2))){
  
    dt_ppv[indx-1+j,i] = ppv[j]
  
    }
}


## keep carrier score ppv values with more than 25% observations in 2000 replicates. 
NA_replicates = apply(dt_ppv, 1, function(x) sum(is.na(x)))
sub_ppv = dt_ppv[NA_replicates < 1500, ]


## calculate means and CI(0.05, 0.95 quantiles) for ppv
means = apply(sub_ppv,1, mean,na.rm = T)
lower_q = apply(sub_ppv, 1, function(x) quantile(x,0.05,na.rm =T))
upper_q = apply(sub_ppv, 1, function(x) quantile(x,0.95,na.rm =T))


x = sequence[NA_replicates < 1500]
s = exp(x)/(1+exp(x))   ## carrier scores on the original scale
y1 = means   ## means of PPV
y2 = lower_q  ## lower bond of ppv 
y3 = upper_q ## upper bond of ppv 


## Code for plot
par(mfrow= c(1,2))

## plot of ppv, with CI
plot(s, y1, type = 'p', pch = 19,cex = 0.1, ylab = "Positive Predictive Value", xlab = "(a) Carrier Scores on Original Scale")
polygon(c(s,rev(s)),c(y2,rev(y3)),col="gray80", border="gray87")
lines(s, y1, type = 'p', pch = 19,cex = 0.1, col = "black")


## plot of ppv, with CI
plot(x, y1, type = 'p', pch = 19,cex = 0.1, ylab = "Positive Predictive Value", xlab = "(b) Carrier Scores on Logit Scale")
polygon(c(x,rev(x)),c(y2,rev(y3)),col="gray80", border="gray87")
lines(x, y1, type = 'p', pch = 19,cex = 0.1, col = "black")
