
## Script to generate Figure3

par(mfrow = c(1,2))

## Vary sample size with theta parameter
setwd("/Users/nicholeyang/Desktop/research2/github code/Data/Vary_SampleSizeTheta")
res1 = read.csv("result_1.csv", header = FALSE)
res2 = read.csv("result_2.csv",header = FALSE)
res3 = read.csv("result_3.csv",header = FALSE)
res4 = read.csv("result_1million_theta.csv",header = FALSE)


boxplot(res1$V1, res2$V1, res3$V1, res4$V1, boxwex=0.4, cex.axis =1, col = "lightgrey", xlab = "(a) Sample Size",ylab =expression(hat(pi)), ylim=c(0,0.3))
mtext("",1,line=1,at=0, cex=0.6)
mtext("1,000",1,line=1,at=1, cex=0.6)
mtext("5,000",1,line=1,at=2, cex=0.6)
mtext("10,000",1,line=1,at=3, cex=0.6)
mtext("1,000,000",1,line=1,at=4, cex=0.6)
abline(h=0.1, col="red", lty= 4)


## Vary theta parameter 
setwd("/Users/nicholeyang/Desktop/research2/github code/Data/Vary_Theta")
theta1 = read.csv("result_theta1.csv", header = FALSE)
theta2 = read.csv("result_theta2.csv", header = FALSE)
theta3 = read.csv("result_theta3.csv", header = FALSE)
theta4 = read.csv("result_theta4.csv", header = FALSE)

boxplot(theta1$V1, theta2$V1, theta3$V1, theta4$V1, boxwex=0.4, cex.axis =1, col = "skyblue", 
        xlab = expression(paste("(b) ", theta)), ylab =expression(hat(pi)),ylim=c(0,0.3))
mtext("",1,line=1,at=0, cex=0.6)
mtext("0.01",1,line=1,at=1, cex=0.6)
mtext("0.05",1,line=1,at=2, cex=0.6)
mtext("0.1",1,line=1,at=3, cex=0.6)
mtext("0.2",1,line=1,at=4, cex=0.6)
abline(h=0.1, col="red", lty= 4)
