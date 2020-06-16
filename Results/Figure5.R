
library(ggplot2)

setwd("/Users/nicholeyang/Desktop/research2/USC_dropbox folder/VUS_paper result")
load("res_summary.RData")


##### BRCApro result
brca = val3[,c("BRCA1or2","GenTestTable.Gene.1.gt","Myriad_panel")]
brca_score =  as.numeric(as.character(brca$BRCA1or2)) 
brca1 = brca[-which(is.na(brca_score)), ]


pos = brca1[which(brca1$Myriad_panel=="Positive"),]
pos_brca = pos$BRCA1or2[which(pos$GenTestTable.Gene.1.gt == "BRCA1" | pos$GenTestTable.Gene.1.gt =="BRCA2")]

VUS = brca1[which(brca1$Myriad_panel=="VUS" | brca1$Myriad_panel=="VUS - Suspect Del/ Likely Pathogenic"),]
VUS_brca = VUS$BRCA1or2[which(VUS$GenTestTable.Gene.1.gt == "BRCA1" | VUS$GenTestTable.Gene.1.gt == "BRCA2")]

neg = brca1[which(brca1$Myriad_panel == "Negative"),]
neg_brca = neg$BRCA1or2


pos1 = as.numeric(as.character(pos_brca))   ###positive carrier probabilities. 76 
VUS1 = as.numeric(as.character(VUS_brca))   #### VUS 30 
neg1 = as.numeric(as.character(neg_brca))    


## to logit scale  
pos2 = log(pos1/(1-pos1))
neg2 = log(neg1/(1-neg1)) 
VUS2 = log(VUS1/(1-VUS1)) 


x = VUS2  #30 
y = pos2  #76
z = neg2  #1066 

sx <- sort(x)
sy <- sort(y)
sz <- sort(z)

sy <- approx(1L:length(y), sy, n = length(x))$y
sz <- approx(1L:length(z), sz, n = length(x))$y



setEPS()
postscript("Figure5.eps")



par(mfrow = c(1,2))

## QQplot
df <- data.frame(sx, sy, sz)
qqplot(df$sx, df$sy, xlim = range(df), ylim = range(df), ylab = "Positive/Negative", pch =20, cex =0.7,xlab = "(a) QQplot for Carrier Score")
mtext("VUS",1,line=2, cex = 0.8)
points(sort(df$sx), sort(df$sz), col = "red", pch =20, cex = 0.7)

## diagonal line
abline(0,1)

legend("topright",legend = c("Positive", "Negative"), col = c(1:2), pch = 20, cex=0.75) 


##
plot(density(neg2, bw = "sj"), col = 2,
     main = "", xlab = "(b) Carrier Score on Logit Scale", xlim = c(-15,10)) ##red
lines(density(pos2, bw = "sj"), col = 1) ##black
lines(density(VUS2, bw = "sj"), col =3)  ##green
legend("topright",legend = c("Positive", "Negative","VUS"), col = c(1:3), pch = 20, cex=0.75)  

dev.off()

