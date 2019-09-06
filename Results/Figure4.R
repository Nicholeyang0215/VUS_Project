
## Script to visualize simulation results of positive predictive value. 

# create sequence of carrier scores, matching the length of row of dt_ppv.
sequence = seq(-9, 10, by = 0.01)

## Load data, dt_ppv
# dt_ppv stores all simulation results of positive predictive value (PPV). 
# Columns of dt_ppv: simulation replicates(1828 in total). 
# Rows of dt_ppv: created based on carrier score of the above sequence. PPVs in each replicate are stored in dt_ppv starting from certain rows according to the carrier score. 
# NAs in a column: in that simulation replicate, we don't have carrier score of that value. 
load("/Users/nicholeyang/Desktop/research2/github code/Data/dt_ppv.RData")


## calculate means of ppv (all replicates)
means = apply(dt_ppv,1, mean,na.rm = T)

## calculate 0.05 and 0.95 quantiles of ppv (all replicates)
lower_q = apply(dt_ppv, 1, function(x) quantile(x,0.05,na.rm =T))
upper_q = apply(dt_ppv, 1, function(x) quantile(x,0.95,na.rm =T))


## remove NAs (Rows that are all NAs)
lower_q2 = lower_q[!is.na(lower_q)]
upper_q2 = upper_q[!is.na(lower_q)]
means2 = means[!is.na(means)]


## get the range of carrier score values with most simulation replicates 
x = sequence[!is.na(means)]

x1 = x[85:1285]
s1 = exp(x[85:1285])/(1+exp(x[85:1285]))   ## carrier scores on the original scale
y1 = means2[85:1285]    ## means of PPV without NAs 
y2 = lower_q2[85:1285]  ## lower bond of PPV without NAs 
y3 = upper_q2[85:1285]  ## upper bond of PPV without NAs



## Code for plot
par(mfrow= c(1,2))

# Plot of PPV on original scale, with CI
plot(s1, y1, type = 'p', pch = 19,cex = 0.1, ylab = "Positive Predictive Value", xlab = "(a) Carrier Probabilities")
polygon(c(s1,rev(s1)),c(y2,rev(y3)),col="gray80", border="gray87")
lines(s1, y1, type = 'p', pch = 19,cex = 0.1, col = "black")


# Plot of PPV on logit scale, with CI
plot(x1, y1, type = 'p', pch = 19,cex = 0.1, ylab = "Positive Predictive Value", xlab = "(b) Logit Carrier Probabilities")
polygon(c(x1,rev(x1)),c(y2,rev(y3)),col="gray80", border="gray87")
lines(x1, y1, type = 'p', pch = 19,cex = 0.1, col = "black")
