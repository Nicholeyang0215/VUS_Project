
## Script to generate Figure 2.


par(mfrow = c(1,3))

## vary sample size 
setwd("/Users/nicholeyang/Desktop/research2/github code/Data/Vary_SampleSize")

res1 = read.csv("result_1.csv", header = FALSE)
res2 = read.csv("result_2.csv",header = FALSE)
res3 = read.csv("result_3.csv",header = FALSE)
res4 = read.csv("result_1million.csv",header = FALSE)


boxplot(res1$V1, res2$V1, res3$V1, res4$V1, boxwex=0.4, cex.axis =1, col = "lightgrey", xlab = "(a) Sample Size", ylab =expression(hat(pi)), ylim=c(0,0.3))
mtext("1000",1,line=1,at=1, cex=0.4)
mtext("5000",1,line=1,at=2, cex=0.4)
mtext("10,000",1,line=1,at=3, cex=0.4)
mtext("1,000,000",1,line=1,at=4, cex=0.4)
abline(h=0.1, col="red", lty= 4)



## vary pi
setwd("/Users/nicholeyang/Desktop/research2/github code/Data/Vary_Pi")
res1 = read.csv("result_1.csv", header = FALSE)
res2 = read.csv("result_2.csv",header = FALSE)
res3 = read.csv("result_3.csv",header = FALSE)
res4 = read.csv("result_4.csv",header = FALSE)

boxplot(res1$V1, res2$V1, res3$V1,res4$V1,  boxwex=0.4, cex.axis =1, col = "skyblue", xlab = expression(paste("(b) True ", pi)),ylab =expression(hat(pi)), ylim=c(0,0.3))
mtext("0.01",1,line=1,at=1, cex=0.4)
mtext("0.05",1,line=1,at=2, cex=0.4)
mtext("0.1",1,line=1,at=3, cex=0.4)
mtext("0.2",1,line=1,at=4, cex=0.4)


## vary P_vus
setwd("/Users/nicholeyang/Desktop/research2/github code/Data/Vary_Pvus")
res1 = read.csv("result_1.csv", header = FALSE)
res2 = read.csv("result_2.csv",header = FALSE)
res3 = read.csv("result_3.csv",header = FALSE)

boxplot(res1$V1, res2$V1, res3$V1, boxwex=0.4, cex.axis =1, col = "khaki1", xlab = expression(paste("(c) ", P(VUS))),ylab =expression(hat(pi)), ylim=c(0,0.3))       
mtext("0.1",1,line=1,at=1, cex=0.4)
mtext("0.2",1,line=1,at=2, cex=0.4)
mtext("0.3",1,line=1,at=3, cex=0.4)
abline(h=0.1, col="red", lty= 4)