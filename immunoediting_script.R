# Author: Clare Weeden

#change the working space and pathway
setwd("/path/to/data")

# read your data
library(readxl)
filename <- "WGS_TRACERx_neoantigen.xlsx"
persample <- read_excel(filename,sheet=6, col_names = TRUE)


library(MASS)
library("writexl")
library("xlsx")

## transform raw TMB (mt/Mb) and # neoantigens by  (SQRT(x+1)-1)*2  ##

#plot xy scatter graph of power_transformed square root data, force 0,0 in robust linear model, 
plot(persample$`PT*sqrt(tmb)`,persample$`PT*sqrt(neoantigens)`,  xlim=c(0,15), ylim=c(0,50),xaxs="i", yaxs="i", bty="L")
ptsqrt_immunoediting.rlm<-(rlm(persample$`PT*sqrt(neoantigens)`~0 + persample$`PT*sqrt(tmb)`))
abline(rlm(persample$`PT*sqrt(neoantigens)`~ 0 + persample$`PT*sqrt(tmb)`))
summary(ptsqrt_immunoediting.rlm)

#data frame of residuals+weights, export as excel
ptsqrt_immunoediting_weights<-data.frame(patientID = persample$Patient, resid = ptsqrt_immunoediting.rlm$residuals, weight = ptsqrt_immunoediting.rlm$w)
write_xlsx(ptsqrt_immunoediting_weights,"Immunoediting_new/ptsqrt_immunoediting_weights.xlsx")

#histo of weights
hist(ptsqrt_immunoediting_weights$weight)
hist(ptsqrt_immunoediting_weights$resid)

#colour as immunoedited if weight less than 1 and negative residual
ptsqrIE<-as.factor(persample$`pt*sqr_immunoedited`)
ptsqrIE.col<-as.factor(persample$`pt*sqr_immunoedited`)
levels(ptsqrIE.col)<-c("red","blue")
plot(persample$`PT*sqrt(tmb)`,persample$`PT*sqrt(neoantigens)`,  xlim=c(0,15), ylim=c(0,50),xaxs="i", yaxs="i", bty="L",xlab = "pt*sqr(TMB)", ylab = "pt*sqr(neoantigen)", pch = 19, col=as.character(ptsqrIE.col))
ptsqrt_immunoediting.rlm<-(rlm(persample$`PT*sqrt(neoantigens)`~0 + persample$`PT*sqrt(tmb)`))
abline(rlm(persample$`PT*sqrt(neoantigens)`~ 0 + persample$`PT*sqrt(tmb)`))
legend("topleft", fill=c("red","blue"), legend = levels(ptsqrIE), cex=0.6)
