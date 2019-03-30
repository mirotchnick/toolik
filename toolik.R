rm(list=ls(all=TRUE))

library(outliers)
library(vegan)
library(MASS)
#library(tree)
library(mgcv)
#library(leaps)
#library(rgl)
#library(Hmisc)
library(picante)
library(nlme)
library(ape)
library(lattice)
#library(lme4)
#library(corrgram)
library(ggplot2)
library(car)
#library(reshape2)
#library(phytools)
#library(Kmisc)
library(arm)
library(broom)

#2008
#setwd("/Users/Nick/Dropbox/Work/Columbia/Spring 2008/Alaska")
setwd("/Users/Nick/Dropbox/Work/Columbia/Spring 2008/Alaska/owen r scripts/redo with break finder")


#toolik <- read.csv("nicksummaryfilteredlikeowen.csv", header = TRUE)
toolik <- read.csv("summaryvariablebreak.csv", header = TRUE)
#toolik3 <- read.csv("master.csv", header = TRUE)


#comparing fixed break point with variable break point data sets

#mean(toolik$Rd.actual)
#mean(toolik2$Rd.actual)

#plot(toolik$Rl/toolik$Rn)
#plot(toolik2$Rd.actual/toolik2$Rn)

#toolik <- aggregate(toolik, list(toolik$Species), mean, na.rm=TRUE)
#toolik2<-aggregate(toolik2,list(toolik2$species),mean,na.rm=TRUE)

#plot(toolik$Rd.actual~toolik2$Rd.actual)




#for (i in 5:ncol(toolik)){
#	toolik[,i] <- as.numeric(as.character(toolik[,i]))
#}





#convert mass-based traits to area-based:

toolik$LMA<- 10000/toolik$SLA
toolik$LMA[313:314]<-mean(toolik$LMA[315:324])
toolik$Rl<-toolik$Rd.actual
toolik$Rl.Rn <- toolik$Rl/toolik$Rn
#toolik$Rl.Rn <- rm.outlier(toolik$Rl.Rn, fill = TRUE)
toolik$Narea <- toolik$N*toolik$LMA
toolik$Parea <- toolik$P*toolik$LMA
toolik$NParea <- toolik$N.P*toolik$LMA
#toolik$Glucose.mg.garea <- toolik$Glucose.mg.g*toolik$LMA
toolik$glucosearea<-toolik$glucose*toolik$LMA
#toolik$Fructose.mg.garea <- toolik$Fructose.mg.g*toolik$LMA
toolik$fructosearea<-toolik$fructose*toolik$LMA
#toolik$Sucrose.mg.garea <- toolik$Sucrose.mg.g*toolik$LMA
toolik$sucrosearea<-toolik$sucrose*toolik$LMA
#toolik$Total.Sugars.mg.garea <- toolik$Total.Sugars.mg.g*toolik$LMA
#toolik$Starch.mg.garea <- toolik$Starch.mg.g*toolik$LMA
toolik$starcharea<-toolik$starch*toolik$LMA
#toolik$TNC.mg.garea <- toolik$TNC.mg.g*toolik$LMA
toolik$TNCarea<-toolik$TNC*toolik$LMA

#convert area-based traits to mass-based:

toolik$Amaxmass<-toolik$Amax/toolik$LMA
toolik$Rlmass <- toolik$Rl/toolik$LMA
toolik$Jsatmass<-toolik$Jsat/toolik$LMA
toolik$Vcsatmass<-toolik$Vcsat/toolik$LMA
toolik$Vosatmass<-toolik$Vosat/toolik$LMA
toolik$Rnmass<-toolik$Rn/toolik$LMA




elev <- toolik[grep("elev",toolik$Replicate),]
toolik <- toolik[grep("amb",toolik$Replicate),]


#attach(toolik)
#attach(elev)


### Data exploration of R-T database ###
# Started on December 13, 2013 #

#data from "R-T Q10 database MH 10.12.2013.xlxs" found in "Database construction" in "Global Q10" on dropbox

#RTdb<-read.delim(pipe("pbpaste"))
#attach(toolik)

#summary(toolik)
#head(toolik)

#### Descriptive statistics and basic exploratory plots ####


#####################################
### Looking for possible outliers ###
#####################################

#ambient outliers:

iqr<-quantile(toolik$Rn, na.rm=T)
toolik$Out<-toolik$Rn<(iqr[2]-2.5*(iqr[4]-iqr[2])) | toolik$Rn > (iqr[4]+2.5*(iqr[4]-iqr[2]))
toolik <- toolik[toolik$Out==FALSE,]

iqr<-quantile(toolik$Rl,na.rm=T)
toolik$Out<-toolik$Rl<(iqr[2]-2.5*(iqr[4]-iqr[2])) | toolik$Rl > (iqr[4]+2.5*(iqr[4]-iqr[2]))
toolik <- toolik[toolik$Out==FALSE,]

iqr<-quantile(toolik$Rl.Rn,na.rm=T)
toolik$Out<-toolik$Rl.Rn<(iqr[2]-2.5*(iqr[4]-iqr[2])) | toolik$Rl.Rn > (iqr[4]+2.5*(iqr[4]-iqr[2]))
toolik <- toolik[toolik$Out==FALSE,]

iqr<-quantile(toolik$Rnmass, na.rm=T)
toolik$Out<-toolik$Rnmass<(iqr[2]-2.5*(iqr[4]-iqr[2])) | toolik$Rnmass > (iqr[4]+2.5*(iqr[4]-iqr[2]))
toolik <- toolik[toolik$Out==FALSE,]

iqr<-quantile(toolik$Rlmass,na.rm=T)
toolik$Out<-toolik$Rlmass<(iqr[2]-2.5*(iqr[4]-iqr[2])) | toolik$Rlmass > (iqr[4]+2.5*(iqr[4]-iqr[2]))
toolik <- toolik[toolik$Out==FALSE,]



#elevated outliers:

iqr<-quantile(elev$Rn, na.rm=T)
elev$Out<-elev$Rn<(iqr[2]-2.5*(iqr[4]-iqr[2])) | elev$Rn > (iqr[4]+2.5*(iqr[4]-iqr[2]))
elev<-elev[elev$Out==FALSE,]

iqr<-quantile(elev$Rl,na.rm=T)
Out<-elev$Rl<(iqr[2]-2.5*(iqr[4]-iqr[2])) | elev$Rl > (iqr[4]+2.5*(iqr[4]-iqr[2]))
elev<-elev[elev$Out==FALSE,]

iqr<-quantile(elev$Rl.Rn,na.rm=T)
elev$Out<-elev$Rl.Rn<(iqr[2]-2.5*(iqr[4]-iqr[2])) | elev$Rl.Rn > (iqr[4]+2.5*(iqr[4]-iqr[2]))
elev<-elev[elev$Out==FALSE,]

iqr<-quantile(elev$Rnmass, na.rm=T)
elev$Out<-elev$Rnmass<(iqr[2]-2.5*(iqr[4]-iqr[2])) | elev$Rnmass > (iqr[4]+2.5*(iqr[4]-iqr[2]))
elev<-elev[elev$Out==FALSE,]

iqr<-quantile(elev$Rlmass,na.rm=T)
Out<-elev$Rlmass<(iqr[2]-2.5*(iqr[4]-iqr[2])) | elev$Rlmass > (iqr[4]+2.5*(iqr[4]-iqr[2]))
elev<-elev[elev$Out==FALSE,]




toolik$QY <- rm.outlier(toolik$QY, fill=TRUE, median=TRUE)
toolik$Cisat <- rm.outlier(toolik$Cisat, fill=TRUE, median=TRUE)







#Rl/Rn outliers above 1:

toolik[Rl.Rn>1,c(1,6,9,54)]



#Borrowed from Keith B "Global respiration data exploration ~ 1st draft.r" as well as the subsequent drafts ##

#Outlier tests will be conducted on respiration values, Q10 calculations, and associated leaf traits #

#Check for normality in variable of interest by plotting histogram#

#First example will use variable "R10a" - area based respiration at 10C #
op<-par(mfrow=c(2,2))
hist(Rn)
qqnorm(Rn)
#histogram appears to be log-normal, will need to transform data: #
#Log transforming data: #

toolik$log.Rn<-log10(toolik$Rn)
hist(toolik$Rn)
qqnorm(toolik$Rn)


op<-par(mfrow=c(2,2))
hist(Rl)
qqnorm(Rl)
#histogram appears to be log-normal, will need to transform data: #
#Log transforming data: #

toolik$log.Rl<-log10(toolik$Rl)
hist(toolik$Rl)
qqnorm(toolik$Rl)


#repeating for mass-based values #
hist(toolik$massRdark)
qqnorm(toolik$massRdark)

#and log-transformed: #
toolik$log.massRdark<-log10(toolik$massRdark)
hist(toolik$log.massRdark)
qqnorm(toolik$log.massRdark)

hist(toolik$massRlight)
qqnorm(toolik$massRlight)

#and log-transformed: #
toolik$log.massRlight<-log10(toolik$massRlight)
hist(toolik$log.massRlight)
qqnorm(toolik$log.massRlight)



# From Keith B: ## Adopting a more objective approach:> one and half times the IQR above the third quartile or below the first quartile

#Graphing quantile cutoffs on histogram

#First up R10a, log transformed:

hist(toolik$Rn, main=paste("Rn (area) ~ All Sites"), sub="Cut offs: 1.5 (red), 2.0 (green), 2.5 (blue)")

#Define quantiles:



#Add cutoff lines on 1.5, 2, and 2.5 time based on above defined quantile of log-transformed data:

abline(v=(iqr[2]-1.5*(iqr[4]-iqr[2])), lty=2, col="red", lwd=2)
abline(v=(iqr[4]+1.5*(iqr[4]-iqr[2])), lty=2, col="red", lwd=2)
abline(v=(iqr[2]-2.0*(iqr[4]-iqr[2])), lty=2, col="green", lwd=2)
abline(v=(iqr[4]+2.0*(iqr[4]-iqr[2])), lty=2, col="green", lwd=2)
abline(v=(iqr[2]-2.5*(iqr[4]-iqr[2])), lty=2, col="blue", lwd=2)
abline(v=(iqr[4]+2.5*(iqr[4]-iqr[2])), lty=2, col="blue", lwd=2)







#again for mass-based log-transformed R10
#remember to redefine quantile:

#hist(toolik$massRdark, main=paste("Rdark (mass) ~ All Sites"), sub="Cut offs: 1.5 (red), 2.0 (green), 2.5 (blue)")
#iqr<-quantile(elev$massRdark, na.rm=T)
#abline(v=(iqr[2]-1.5*(iqr[4]-iqr[2])), lty=2, col="red", lwd=2)
#abline(v=(iqr[4]+1.5*(iqr[4]-iqr[2])), lty=2, col="red", lwd=2)
#abline(v=(iqr[2]-2.0*(iqr[4]-iqr[2])), lty=2, col="green", lwd=2)
#abline(v=(iqr[4]+2.0*(iqr[4]-iqr[2])), lty=2, col="green", lwd=2)
#abline(v=(iqr[2]-2.5*(iqr[4]-iqr[2])), lty=2, col="blue", lwd=2)
#abline(v=(iqr[4]+2.5*(iqr[4]-iqr[2])), lty=2, col="blue", lwd=2)

#Save plot in "Outliers" file on DropBox

#How many outliers are 2.5x (or 2.0, or 1.5x) beyond keeping?#
##Using Keith B's definition of outliers (here with log-transformed R10-mass based at a 2.5 level cutoff:##

#Out<-elev$massRdark<(iqr[2]-2.5*(iqr[4]-iqr[2])) | elev$massRdark > (iqr[4]+2.5*(iqr[4]-iqr[2]))

## How many outliers? ##
length(which(Out))

## Which values are they? Gives specific row numbers of the outlier values. ##

which(Out)

## Currently plan to manually remove individual outlier value from Excel spreadsheet, not whole row.#
## Need to have a w and w/o outlier spreadsheet for keeping ##

elev <- cbind(elev,Out)
elev <- elev[elev$Out==FALSE,]
toolik$Out <- Out
toolik <- toolik[toolik$Out==FALSE,]


elev$Amax <- rm.outlier(elev$Amax, fill=TRUE, median=TRUE)
elev$QY <- rm.outlier(elev$QY, fill=TRUE, median=TRUE)
elev$QY <- rm.outlier(elev$QY, fill=TRUE, median=TRUE)
elev$QY <- rm.outlier(elev$QY, fill=TRUE, median=TRUE)
elev$QY <- rm.outlier(elev$QY, fill=TRUE, median=TRUE)
elev$QY <- rm.outlier(elev$QY, fill=TRUE, median=TRUE)

#repeat as necessary:
a <- rm.outlier(elev$Rl.Rn, fill=TRUE, median=TRUE)
a <- rm.outlier(a, fill=TRUE, median=TRUE)
mean(a)
plot(a)
elev$Rl.Rn <- a


attach(elev)

#elev <- aggregate(elev, list(elev$Species), mean, na.rm=TRUE)
#attach(elev)

b <- rm.outlier(toolik$Rl.Rn, fill=TRUE, median=TRUE)



# write.csv(toolik, "nicksummaryedit.csv")
# write.csv(amb, "nicksummaryeditamb.csv")
# write.csv(elev, "nicksummaryeditelev.csv")











#display raw data


# hist(1-Rl.Rn)
attach(toolik)







ggplot(as.data.frame(species), aes(Rl.Rn, fill = species)) + geom_density(alpha = 0.3) + theme(legend.position="top") + labs(x="Rl/Rn")

ggplot(as.data.frame(Functional), aes(Rl.Rn, fill = Functional)) + geom_density(alpha = 0.3) + labs(x="Rl/Rn", fill="Functional group")


plot(Rl.Rn, xlab="individuals", ylab=expression(R[L]/R[D]~"( %)"))
abline(a=mean(Rl.Rn), b=0,)
abline(a=mean(Rl.Rn)+(sd(Rl.Rn)/sqrt(length(Rl.Rn)))*qnorm(0.975), b=0, lty="dotted")
abline(a=mean(Rl.Rn)-(sd(Rl.Rn)/sqrt(length(Rl.Rn)))*qnorm(0.975), b=0, lty="dotted")


plot(Rl.Rn,col=as.numeric(species),pch=as.numeric(species), xlab="individuals", ylab=expression(R[L]/R[D]~"( %)"))
abline(a=mean(Rlight.Rdark400), b=0,)
abline(a=mean(Rlight.Rdark400)+(sd(Rlight.Rdark400)/sqrt(length(Rlight.Rdark400)))*qnorm(0.975), b=0, lty="dotted")
abline(a=mean(Rlight.Rdark400)-(sd(Rlight.Rdark400)/sqrt(length(Rlight.Rdark400)))*qnorm(0.975), b=0, lty="dotted")




means <- aggregate(toolik, list(toolik$species), mean, na.rm=TRUE)

ci <- sd(toolik$Rl.Rn)/sqrt(length(toolik$Rl.Rn))*qnorm(0.975)

ggplot(as.data.frame(toolik), aes(reorder(species,Rl.Rn,function(x)mean(x)), Rl.Rn)) + stat_summary() +theme(text=element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) + theme(plot.margin=unit(c(1,1,1,2),"cm")) + labs(x="", y="Rl/Rd (%)") + theme(panel.background = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black")) + geom_hline(yintercept=mean(toolik$Rl.Rn), color="blue") + geom_hline(yintercept=mean(toolik$Rl.Rn)+ci, color="grey", linetype=2) + geom_hline(yintercept=mean(toolik$Rl.Rn)-ci, color="grey", linetype=2)


ggplot(as.data.frame(toolik), aes(reorder(species,Rl.Rn,function(x)mean(x)), Rl.Rn)) + stat_summary(aes(y=-Rn,group=1,color="x",labs="Rn"), geom="line") + scale_color_manual(values=c("x"="magenta","y"="cyan"),labels=c("Rd","Rl")) + stat_summary(aes(y=-Rn,group=1), geom="point") + stat_summary(aes(y=-Rn), fun.data = mean_se, geom = "pointrange") + stat_summary(aes(y=-Rl,group=1,color="y"),geom="line") + labs(color="") + stat_summary(aes(y=-Rl,group=1), geom="point") + stat_summary(aes(y=-Rl), fun.data = mean_se, geom = "pointrange") + theme(text=element_text(size=12), axis.text.x = element_text(angle = 45, hjust = 1, size=12)) + theme(plot.margin=unit(c(1,1,1,2),"cm")) + labs(x="", y="R") + theme(legend.position="top") + theme(panel.background = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"))

#plot R variance by species
boxplot(Rl.Rn~Species, las=2, ylab=expression(R[L]/R[D]~"( %)"), par(oma=c(5,1,1,1)), cex.axis=0.75)

boxplot(Rl~Species, las=2, ylab=expression(R[L]), par(oma=c(5,1,1,1)), cex.axis=0.75)

boxplot(Rn~Species, las=2, ylab=expression(R[D]), par(oma=c(5,1,1,1)), cex.axis=0.75)


#toolik$Rl <- rm.outlier(toolik$Rl, fill=TRUE, median=TRUE)
#toolik$Rn <- rm.outlier(toolik$Rn, fill=TRUE, median=TRUE)
#attach(toolik)

plot(Rl~Rn, pch=16, cex=0.5, xlim=c(-3,-0.4), ylim=c(-3,-0.4))
abline(a=0,b=1)
text(Rn, Rl, labels=Group.1, cex=0.4, adj=-0.1)

#plots for elevated CO2
elev <- aggregate(elev, list(Species), mean, na.rm=TRUE)
attach(elev)
plot(Rl~Rn, pch=16, cex=0.5)
abline(a=0,b=1)
text(Rn, Rl, labels=Group.1, cex=0.4, adj=-0.1)





plot(Rlight.Rdark400,col=as.numeric(Functional),pch=as.numeric(Functional), xlab="individuals", ylab=expression(R[L]/R[D]~"( %)"))
legend(-5,41, levels(Functional), pch=1:5, col=1:5)
abline(a=mean(Rlight.Rdark400), b=0,)
abline(a=mean(Rlight.Rdark400)+(sd(Rlight.Rdark400)/sqrt(length(Rlight.Rdark400)))*qnorm(0.975), b=0, lty="dotted")
abline(a=mean(Rlight.Rdark400)-(sd(Rlight.Rdark400)/sqrt(length(Rlight.Rdark400)))*qnorm(0.975), b=0, lty="dotted")

#anovas of functional groups

#area
an1 <- oneway.test(toolik$Rl.Rn~as.factor(toolik$Functional))
an1
an2 <- oneway.test(toolik$Rl.Rn~as.factor(toolik$PFTLPJ))
an2
an3 <- oneway.test(toolik$Rl.Rn~as.factor(toolik$PFTSheffield))
an3
an4 <- oneway.test(toolik$Rl.Rn~as.factor(toolik$PFTJULES))
an4

an5 <- oneway.test(toolik$Rl~as.factor(toolik$Functional))
an5
an6 <- oneway.test(toolik$Rl~as.factor(toolik$PFTLPJ))
an6
an7 <- oneway.test(toolik$Rl~as.factor(toolik$PFTSheffield))
an7
an8 <- oneway.test(toolik$Rl~as.factor(toolik$PFTJULES))
an8

an9 <- oneway.test(toolik$Rn~as.factor(toolik$Functional))
an9
an10 <- oneway.test(toolik$Rn~as.factor(toolik$PFTLPJ))
an10
an11 <- oneway.test(toolik$Rn~as.factor(toolik$PFTSheffield))
an11
an12 <- oneway.test(toolik$Rn~as.factor(toolik$PFTJULES))
an12


#mass
an1 <- oneway.test(toolik$Rl.Rn~as.factor(toolik$Functional))
an1
an2 <- oneway.test(toolik$Rl.Rn~as.factor(toolik$PFTLPJ))
an2
an3 <- oneway.test(toolik$Rl.Rn~as.factor(toolik$PFTSheffield))
an3
an4 <- oneway.test(toolik$Rl.Rn~as.factor(toolik$PFTJULES))
an4

an5 <- oneway.test(toolik$Rlmass~as.factor(toolik$Functional))
an5
# an6 <- oneway.test(toolik$Rlmass~as.factor(toolik$PFTLPJ))
an6
an7 <- oneway.test(toolik$Rlmass~as.factor(toolik$PFTSheffield))
an7
an8 <- oneway.test(toolik$Rlmass~as.factor(toolik$PFTJULES))
an8

an9 <- oneway.test(toolik$Rnmass~as.factor(toolik$Functional))
an9
an10 <- oneway.test(toolik$Rnmass~as.factor(toolik$PFTLPJ))
an10
an11 <- oneway.test(toolik$Rnmass~as.factor(toolik$PFTSheffield))
an11
an12 <- oneway.test(toolik$Rnmass~as.factor(toolik$PFTJULES))
an12

boxplot(Rl.Rn~Functional, las=2, ylab=expression(R[L]/R[D]~"( %)"), par(oma=c(3,1,1,1)))

boxplot(Rlmass~Functional, las=2, ylab=expression(R[L]), par(oma=c(3,1,1,1)))

boxplot(Rnmass~Functional, las=2, ylab=expression(R[D]), par(oma=c(3,1,1,1)))

boxplot(Rl.Rn~PFTLPJ, las=2, ylab=expression(R[L]/R[D]~"( %)"), par(oma=c(3,1,1,1)))

boxplot(Rlmass~PFTLPJ, las=2, ylab=expression(R[L]/R[D]~"( %)"), par(oma=c(3,1,1,1)))

boxplot(Rnmass~PFTLPJ, las=2, ylab=expression(R[L]/R[D]~"( %)"), par(oma=c(3,1,1,1)))

boxplot(Rl.Rn~ PFTSheffield, las=2, ylab=expression(R[L]/R[D]~"( %)"), par(oma=c(3,1,1,1)))

boxplot(Rlmass~ PFTSheffield, las=2, ylab=expression(R[L]/R[D]~"( %)"), par(oma=c(3,1,1,1)))

boxplot(Rnmass~ PFTSheffield, las=2, ylab=expression(R[L]/R[D]~"( %)"), par(oma=c(3,1,1,1)))

boxplot(Rl.Rn~ PFTJULES, las=2, ylab=expression(R[L]/R[D]~"( %)"), par(oma=c(3,1,1,1)))

boxplot(Rlmass~ PFTJULES, las=2, ylab=expression(R[L]/R[D]~"( %)"), par(oma=c(3,1,1,1)))

boxplot(Rnmass~ PFTJULES, las=2, ylab=expression(R[L]/R[D]~"( %)"), par(oma=c(3,1,1,1)))



plotCI(mean(Rlight.Rdark400))

plot(Rlight.Rdark400,col=as.numeric(Family),pch=as.numeric(Family), xlab="individuals", ylab=expression(R[L]/R[D]~"( %)"))
legend(-8,62, levels(Family), pch=1:5, col=1:5, bty="n", cex=0.9)
abline(a=mean(Rlight.Rdark400), b=0,)
abline(a=mean(Rlight.Rdark400)+(sd(Rlight.Rdark400)/sqrt(length(Rlight.Rdark400)))*qnorm(0.975), b=0, lty="dotted")
abline(a=mean(Rlight.Rdark400)-(sd(Rlight.Rdark400)/sqrt(length(Rlight.Rdark400)))*qnorm(0.975), b=0, lty="dotted")
anova <- oneway.test(toolik$Rlight.Rdark400~as.factor(toolik$Family))
boxplot(Rlight.Rdark400~Family, las=2, ylab=expression(R[L]/R[D]~"( %)"), par(oma=c(2,1,1,1)))


#traits

alltraits <- as.data.frame(cbind(Species,avgTleaf,Amax,QY,PSsat,Jsat,Vcsat,Vosat,Gssat,Cisat,PS100,J100,Vc100,Vo100,Gs100,Ci100,Rn,Rl,Rl.Rn,LI_ratio,Rlight.Asat,Rdark.Asat,massAmax,massAsat,massJsat,massVcsat,massVosat,massJ100,massVc100,massVo100,massRlight,massRdark,NAmax,NAsat,NRdark,PAmax,PAsat,PRdark,N,P,NP,Narea,Parea,NParea,LMA,Glucose.mg.garea,Fructose.mg.garea,Sucrose.mg.garea,Total.Sugars.mg.garea,Starch.mg.garea,TNC.mg.garea,Glucose.mg.g,Fructose.mg.g,Sucrose.mg.g,Total.Sugars.mg.g,Starch.mg.g,TNC.mg.g))

alltraits <- as.data.frame(cbind(Amax,QY,PSsat,Jsat,Vcsat,Vosat,Gssat,Cisat,Narea,Parea,LMA,Amaxmass,Jsatmass,Vcsatmass,Vosatmass,N,P,Rl,Rn,Rlmass,Rnmass,Rl.Rn))
alltraits <- decostand(alltraits, "standardize")
attach(alltraits)

#areatraits <- as.data.frame(cbind(avgTleaf,Amax,QY,PSsat,Jsat,Vcsat,Vosat,Gssat,Cisat,PS100,J100,Vc100,Vo100,Gs100,Ci100,Rn,LI_ratio,Rlight.Asat,Rdark.Asat,NAmax,NAsat,NRdark,PAmax,PAsat,PRdark,Narea,Parea,NParea,LMA,Glucose.mg.garea,Fructose.mg.garea,Sucrose.mg.garea,Total.Sugars.mg.garea,Starch.mg.garea,TNC.mg.garea))

#areanonresp <- as.data.frame(cbind(avgTleaf,Amax,QY,PSsat,Jsat,Vcsat,Vosat,Gssat,Cisat,PS100,J100,Vc100,Vo100,Gs100,Ci100,LI_ratio,NAmax,NAsat,PAmax,PAsat,Narea,Parea,NParea,LMA,Glucose.mg.garea,Fructose.mg.garea,Sucrose.mg.garea,Total.Sugars.mg.garea,Starch.mg.garea,TNC.mg.garea))

#stripped down area traits
#areanonresp <-as.data.frame(cbind(Amax,QY,PSsat,Jsat,Vcsat,Vosat,Gssat,Cisat,PS100,J100,Vc100,Vo100,Gs100,Ci100,LI_ratio,Narea,Parea,NParea,LMA,Glucose.mg.garea,Fructose.mg.garea,Sucrose.mg.garea,Total.Sugars.mg.garea,Starch.mg.garea,TNC.mg.garea))


#masstraits <- cbind(avgTleaf,Amax,QY,PSsat,Jsat,Vcsat,Vosat,Gssat,Cisat,PS100,J100,Vc100,Vo100,Gs100,Ci100,Rn,LI_ratio,Rlight.Asat,Rdark.Asat,massAmax,massAsat,massJsat,massVcsat,massVosat,massJ100,massVc100,massVo100,massRlight,massRdark,N,P,NP)

#massnonresp <- as.data.frame(cbind(avgTleaf,QY,PSsat,Jsat,Vcsat,Vosat,Gssat,Cisat,PS100,J100,Vc100,Vo100,Gs100,Ci100,LI_ratio,massAmax,massAsat,massJsat,massVcsat,massVosat,massJ100,massVc100,massVo100,N,P,NP))

#stripped down mass traits
#massnonresp <- as.data.frame(cbind(QY,PSsat,Gssat,Cisat,PS100,Gs100,Ci100,LI_ratio,massAmax,massAsat,massJsat,massVcsat,massVosat,massJ100,massVc100,massVo100,N,P,NP,Glucose.mg.g,Fructose.mg.g,Sucrose.mg.g,Total.Sugars.mg.g,Starch.mg.g,TNC.mg.g))

rltraits <- cbind(Rl,Rlight.Asat,massRlight,NRlight,PRlight,Rlmass)

rdarktraits <- cbind(Rd.apparent,Rd.actual,Rn,massRdark)




toolik[,13:length(toolik)] <- decostand(toolik[,13:length(toolik)-1], "standardize")
attach(toolik)

elev[,13:length(elev)-1] <- decostand(elev[,13:length(elev)-1], "standardize")
attach(elev)

areanonresp <-as.data.frame(cbind(Amax,QY,PSsat,Jsat,Vcsat,Vosat,Gssat,Cisat,Narea,Parea,LMA,glucosearea,fructosearea,sucrosearea,starcharea))

massnonresp <- as.data.frame(cbind(Amaxmass,QY,PSsat,Jsatmass,Vcsatmass,Vosatmass,Gssat,Cisat,N,P,LMA,glucose,fructose,sucrose,starch))






toolik <- aggregate(toolik, list(species), mean, na.rm=TRUE)
attach(toolik)





attach(elev)

toolik400 <- (toolik[,5:35])
toolik400 <- (toolik[,6:12])

rtable <- data.frame(as.character(Group.1),Rl,Rn,Rl.Rn)
rtable
write.csv(rtable,"rtable.csv")

#standardize data
toolikstand <- decostand(toolik400, "range")
toolikstand <- decostand(toolik[,9:length(toolik)], "standardize")
attach(toolikstand)

#traits <- toolikstand[,1:6]
respiration <- toolik400[7:10]
#or
respiration <- toolik400[9:11]

plot(SLA ~ Species)

cor(toolik400)

pca <- princomp(alltraits, na.action=na.omit)
summary(pca,loadings=TRUE,cutoff=.001)
plot(pca)
biplot(pca)

pca2 <- princomp(areatraits, na.action=na.omit)
summary(pca2, loadings=TRUE, cutoff=0.001)
plot(pca2)
biplot(pca2)

pca3 <- princomp(masstraits, na.action=na.omit)
summary(pca3, loadings=TRUE, cutoff=0.001)
plot(pca3)
biplot(pca3)

#pca4 performs a principal components analysis of all area-based traits that are not rates of respiration

pca4 <- princomp(areanonresp, na.action=na.omit)
summary(pca4, loadings=TRUE, cutoff=0.001)
plot(pca4)
biplot(pca4,cex=0.75)
ggord(pca4)
biplot(pca4, choices=2:3)
pca4scores <- as.data.frame(pca4$scores)
pcamodel4 <- lm(Rl~pca4scores$Comp.1+pca4scores$Comp.2+pca4scores$Comp.3+pca4scores$Comp.4)
summary(pcamodel4)
summary(step(pcamodel4))
pov <- pca4$sdev^2/sum(pca4$sdev^2)
cpov <- cumsum(pov)
var <- rbind("proportion of variance"=pov, "cumulative proportion" = cpov)
write.csv(round(var, digits=4), "pca4var.csv")

traitnames <- as.vector(names(pca4$loadings[,1]))
x <- barplot(abs(pca4$loadings[,1]), xaxt="n")
text(x, par("usr")[3], labels=traitnames, xpd=TRUE, srt=45, adj=c(1.1,1.1), cex=0.6)
load <- cbind(names(pca4$loadings[,1]),pca4$loadings[,1])
sort <- load[order(abs(as.numeric(load[,2])), decreasing=TRUE),]
sort
data <- names(sort[,2])
clip <- pipe("pbcopy", "w")
write.table(sort, file=clip, quote=FALSE, row.names=FALSE)
close(clip)
table<-sort
lm(Rl~pca4scores$Comp.1)
summary(lm(Rl~pca4scores$Comp.1))
plot(Rl~pca4scores$Comp.1, xlab="principal component 1", ylab=expression(paste(R[L])))
abline(lm(Rl~pca4scores$Comp.1))
rsq <- round(summary(lm(Rl~pca4scores$Comp.1))$adj.r.squared, digits=2)
text(-9.5,1.3,expression(R^2==rsq))


traitnames <- as.vector(names(pca4$loadings[,1]))
x <- barplot(abs(pca4$loadings[,2]), xaxt="n")
text(x, par("usr")[3], labels=traitnames, xpd=TRUE, srt=45, adj=c(1.1,1.1), cex=0.6)
load <- cbind(names(pca4$loadings[,1]),pca4$loadings[,2])
sort <- load[order(abs(as.numeric(load[,2])), decreasing=TRUE),]
sort
data <- names(sort[,2])
clip <- pipe("pbcopy", "w")
write.table(data, file=clip, quote=FALSE, row.names=FALSE)
close(clip)
table<-cbind(table,sort)
summary(lm(Rl~pca4scores$Comp.2))


traitnames <- as.vector(names(pca4$loadings[,3]))
x <- barplot(abs(pca4$loadings[,3]), xaxt="n")
text(x, par("usr")[3], labels=traitnames, xpd=TRUE, srt=45, adj=c(1.1,1.1), cex=0.6)
load <- cbind(names(pca4$loadings[,1]),pca4$loadings[,3])
sort <- load[order(abs(as.numeric(load[,2])), decreasing=TRUE),]
sort
data <- names(sort[,2])
clip <- pipe("pbcopy", "w")
write.table(data, file=clip, quote=FALSE, row.names=FALSE)
close(clip)
table<-cbind(table,sort)
summary(lm(Rl~pca4scores$Comp.3))
plot(Rl~pca4scores$Comp.3, xlab="principal component 3", ylab=expression(paste(R[L])))
abline(lm(Rl~pca4scores$Comp.3))
rsq <- round(summary(lm(Rl~pca4scores$Comp.3))$adj.r.squared, digits=2)
text(-9.5,1.3,expression(R^2==rsq))



traitnames <- as.vector(names(pca4$loadings[,4]))
x <- barplot(abs(pca4$loadings[,4]), xaxt="n")
text(x, par("usr")[3], labels=traitnames, xpd=TRUE, srt=45, adj=c(1.1,1.1), cex=0.6)
load <- cbind(names(pca4$loadings[,1]),pca4$loadings[,4])
sort <- load[order(abs(as.numeric(load[,2])), decreasing=TRUE),]
sort
data <- names(sort[,2])
clip <- pipe("pbcopy", "w")
write.table(data, file=clip, quote=FALSE, row.names=FALSE)
close(clip)
table<-cbind(table,sort)
write.csv(table,"pca4loadings.csv")
summary(lm(Rl~pca4scores$Comp.4))

summary(lm(Rl~Jsat))
plot(Rl~Jsat)

summary(lm(Rl~Vcsat))
plot(Rl~Vcsat)	

summary(lm(Rl~Amax))
plot(Rl~Amax)
abline(lm(Amax~Rl))

summary(lm(Rl~Narea))
plot(Rl~Narea)

summary(lm(Rl~Gssat))
plot(Rl~Gssat)

summary(lm(Rl~Gs100))
plot(Rl~Gs100)

summary(lm(Rl~Vo100))
plot(Rl~Vo100)

summary(lm(Rl~Cisat))
plot(Rl~Cisat)

summary(lm(Rl~J100))
plot(Rl~J100)

summary(lm(Rl~Vc100))
plot(Rl~Vc100)

summary(lm(Rl~LMA))

summary(lm(Rl~Total.Sugars.mg.garea))

summary(lm(Rl~TNC.mg.garea))

summary(lm(Rl~Parea))

summary(lm(Rl~starcharea))
plot(Rl~starcharea)

summary(lm(Rl~sucrosearea))

summary(lm(Rl~QY))

summary(lm(Rl~PSsat))









#pca5 performs a principal components analysis of all mass-based traits that are not rates of respiration

pca5 <- princomp(massnonresp, na.action=na.omit)
summary(pca5, loadings=TRUE, cutoff=0.001)
plot(pca5)
biplot(pca5, cex=0.75)
biplot(pca5, choices=2:3)
pca5scores <- as.data.frame(pca5$scores)
pcamodel5 <- lm(Rlmass~pca5scores$Comp.1+pca5scores$Comp.2+pca5scores$Comp.3+pca5scores$Comp.4)
summary(pcamodel5)
summary(step(pcamodel5))
pov <- pca5$sdev^2/sum(pca5$sdev^2)
cpov <- cumsum(pov)
var <- rbind("proportion of variance"=pov, "cumulative proportion" = cpov)
write.csv(round(var, digits=4), "pca5var.csv")

traitnames <- as.vector(names(pca5$loadings[,1]))
x <- barplot(abs(pca5$loadings[,1]), xaxt="n")
text(x, par("usr")[3], labels=traitnames, xpd=TRUE, srt=45, adj=c(1.1,1.1), cex=0.6)
load <- cbind(names(pca5$loadings[,1]),pca5$loadings[,1])
sort <- load[order(abs(as.numeric(load[,2])), decreasing=TRUE),]
sort
data <- names(sort[,2])
clip <- pipe("pbcopy", "w")
write.table(data, file=clip, quote=FALSE, row.names=FALSE)
close(clip)
table<-sort
summary(lm(Rlmass~pca5scores$Comp.1))

traitnames <- as.vector(names(pca5$loadings[,1]))
x <- barplot(abs(pca5$loadings[,2]), xaxt="n")
text(x, par("usr")[3], labels=traitnames, xpd=TRUE, srt=45, adj=c(1.1,1.1), cex=0.6)
load <- cbind(names(pca5$loadings[,1]),pca5$loadings[,2])
sort <- load[order(abs(as.numeric(load[,2])), decreasing=TRUE),]
sort
data <- names(sort[,2])
clip <- pipe("pbcopy", "w")
write.table(data, file=clip, quote=FALSE, row.names=FALSE)
close(clip)
table<-cbind(table,sort)
summary(lm(Rlmass~pca5scores$Comp.2))

traitnames <- as.vector(names(pca5$loadings[,1]))
x <- barplot(abs(pca5$loadings[,3]), xaxt="n")
text(x, par("usr")[3], labels=traitnames, xpd=TRUE, srt=45, adj=c(1.1,1.1), cex=0.6)
load <- cbind(names(pca5$loadings[,1]),pca5$loadings[,3])
sort <- load[order(abs(as.numeric(load[,2])), decreasing=TRUE),]
sort
data <- names(sort[,2])
clip <- pipe("pbcopy", "w")
write.table(data, file=clip, quote=FALSE, row.names=FALSE)
close(clip)
table<-cbind(table,sort)
summary(lm(Rlmass~pca5scores$Comp.3))

traitnames <- as.vector(names(pca5$loadings[,1]))
x <- barplot(abs(pca5$loadings[,4]), xaxt="n")
text(x, par("usr")[3], labels=traitnames, xpd=TRUE, srt=45, adj=c(1.1,1.1), cex=0.6)
load <- cbind(names(pca5$loadings[,1]),pca5$loadings[,4])
sort <- load[order(abs(as.numeric(load[,2])), decreasing=TRUE),]
sort
data <- names(sort[,2])
clip <- pipe("pbcopy", "w")
write.table(data, file=clip, quote=FALSE, row.names=FALSE)
close(clip)
table<-cbind(table,sort)
write.csv(table,"pca5loadings.csv")
summary(lm(Rlmass~pca5scores$Comp.4))

summary(lm(Rlmass~Vcsatmass))
plot(Rlmass ~massVcsat)

summary(lm(Rlmass ~Jsatmass))
plot(Rlmass ~massJsat)

summary(lm(Rlmass ~Amax))

summary(lm(Rlmass ~massAsat))
plot(Rlmass ~massAsat)

summary(lm(Rlmass ~N))
plot(Rlmass ~N)

summary(lm(Rlmass ~Vo100))

summary(lm(massRlight~J100))

summary(lm(massRlight~Vc100))

summary(lm(Rlmass ~LMA))

summary(lm(Rlmass ~Total.Sugars.mg.g))
plot(Rlmass ~Total.Sugars.mg.g)

summary(lm(Rlmass ~Glucose.mg.g))
plot(Rlmass ~Glucose.mg.g)

summary(lm(Rlmass ~P))
plot(Rlmass ~P)

summary(lm(Rlmass ~NP))

summary(lm(Rlmass ~starch))
plot(Rlmass~starch)

summary(lm(Rlmass ~sucrose))
plot(Rlmass ~Sucrose.mg.g)

summary(lm(Rlmass ~Fructose.mg.g))
plot(Rlmass ~Fructose.mg.g)

summary(lm(massRlight~TNC.mg.g))
plot(massRlight~TNC.mg.g)

summary(lm(massRlight~Gssat))
plot(massRlight~Gssat)

summary(lm(massRlight~Gs100))
plot(massRlight~Gs100)

summary(lm(Rlmass~QY))
plot(massRlight~QY)

summary(lm(Rlmass~PSsat))

summary(lm(massRlight~Vcsat))

summary(lm(massRlight~Vosat))
plot(massRlight~Vosat)







pcamodel6 <- lm(Rl.Rn~pca4scores$Comp.1+pca4scores$Comp.2+pca4scores$Comp.3+pca4scores$Comp.4)
summary(pcamodel6)
summary(step(pcamodel6))

summary(lm(Rl.Rn~pca4scores$Comp.1))
summary(lm(Rl.Rn~pca4scores$Comp.2))
summary(lm(Rl.Rn~pca4scores$Comp.3))
summary(lm(Rl.Rn~pca4scores$Comp.4))

summary(lm(Rl.Rn~sucrose))
plot(Rl.Rn~Sucrose.mg.g)

summary(lm(Rl.Rn~TNC.mg.garea))
plot(Rl.Rn~TNC.mg.g)

summary(lm(Rl.Rn~LMA))
plot(Rl.Rn~LMA)

summary(lm(Rl.Rn~Jsat))

summary(lm(Rl.Rn~Vcsat))

summary(lm(Rl.Rn~Amax))

summary(lm(Rl.Rn~Narea))

summary(lm(Rl.Rn~Vo100))

summary(lm(Rl.Rn~J100))

summary(lm(Rl.Rn~Vc100))

summary(lm(Rl.Rn~Total.Sugars.mg.garea))

summary(lm(Rl.Rn~Parea))

summary(lm(Rl.Rn~starch))

summary(lm(Rl.Rn~QY))

summary(lm(Rl.Rn~PSsat))





pcamodel8 <- lm(Rnmass~pca5scores$Comp.1+pca5scores$Comp.2+pca5scores$Comp.3+pca5scores$Comp.4)
summary(pcamodel8)
summary(step(pcamodel8))

summary(lm(Rnmass ~pca5scores$Comp.1))
summary(lm(Rnmass ~pca5scores$Comp.2))
summary(lm(Rnmass ~pca5scores$Comp.3))
summary(lm(Rnmass ~pca5scores$Comp.4))

summary(lm(massRdark~TNC.mg.g))
plot(massRdark~TNC.mg.g)

summary(lm(massRdark~Total.Sugars.mg.g))
plot(massRdark~Total.Sugars.mg.g)

summary(lm(Rnmass~N))
plot(massRdark~N)

summary(lm(Rnmass~Jsat))
plot(massRdark~Jsat)

summary(lm(Rnmass~Vcsat))
plot(massRdark~Vcsat)

summary(lm(Rnmass~Gssat))
plot(massRdark~Gssat)

summary(lm(Rnmass~Amax))

summary(lm(massRdark~Vc100))

summary(lm(massRdark~Vo100))

summary(lm(massRdark~J100))

summary(lm(Rnmass~LMA))

summary(lm(Rnmass~P))

summary(lm(Rnmass~starch))

summary(lm(Rnmass~sucrose))

summary(lm(Rnmass~QY))

summary(lm(Rnmass~PSsat))



pcamodel7 <- lm(Rl.Rn~pca5scores$Comp.1+pca5scores$Comp.2+pca5scores$Comp.3+pca5scores$Comp.4)
summary(pcamodel7)
summary(step(pcamodel7))

summary(lm(Rl.Rn~pca5scores$Comp.1))
summary(lm(Rl.Rn~pca5scores$Comp.2))
summary(lm(Rl.Rn~pca5scores$Comp.3))
summary(lm(Rl.Rn~pca5scores$Comp.4))

summary(lm(Rl.Rn~N))

summary(lm(Rl.Rn~TNC.mg.g))

summary(lm(Rl.Rn~P))

summary(lm(Rl.Rn~Starch.mg.g))

summary(lm(Rl.Rn~Sucrose.mg.g))

summary(lm(Rl.Rn~Glucose.mg.g))

summary(lm(Rl.Rn~Fructose.mg.g))

summary(lm(Rl.Rn~QY))

summary(lm(Rl.Rn~Jsat))

summary(lm(Rl.Rn~J100))

summary(lm(Rl.Rn~Vcsat))

summary(lm(Rl.Rn~Vc100))

summary(lm(Rl.Rn~Amax))

summary(lm(Rl.Rn~Vo100))

summary(lm(Rl.Rn~LMA))

summary(lm(Rl.Rn~Total.Sugars.mg.g))

summary(lm(Rl.Rn~PSsat))




pcamodel9 <- lm(Rn~pca4scores$Comp.1+pca4scores$Comp.2+pca4scores$Comp.3+pca4scores$Comp.4)
summary(pcamodel9)
summary(step(pcamodel9))

summary(lm(Rn~pca4scores$Comp.1))
summary(lm(Rn~pca4scores$Comp.2))
summary(lm(Rn~pca4scores$Comp.3))
summary(lm(Rn~pca4scores$Comp.4))

summary(lm(Rn~LMA))
plot(Rn~LMA)

summary(lm(Rn~Total.Sugars.mg.garea))
plot(Rn~Total.Sugars.mg.garea)

summary(lm(Rn~TNC.mg.garea))
plot(Rn~TNC.mg.garea)

summary(lm(Rn~Gs100))
plot(Rn~Gs100)

summary(lm(Rn~Gssat))
plot(Rn~Gssat)

summary(lm(Rn~Vo100))
plot(Rn~Vo100)

summary(lm(Rn~Vcsat))
plot(Rn~Vcsat)

summary(lm(Rn~Vc100))

summary(lm(Rn~Jsat))
plot(Rn~Jsat)

summary(lm(Rn~J100))

summary(lm(toolik$Rn~toolik$PSsat))
plot(toolik$Rn~toolik$PSsat)

summary(lm(Rn~Amax))
plot(Rn~Amax)

summary(lm(Rn~starcharea))
plot(Rn~Starch.mg.garea)

summary(lm(Rn~Parea))
plot(Rn~Parea)

summary(lm(Rn~Narea))
plot(Rn~Narea)

summary(lm(Rl~sucrosearea))

summary(lm(Rl~QY))

summary(lm(Rl~PSsat))





pcamodel10 <- lm(Rl~Amax+QY+PSsat+Jsat+Vcsat+Vosat+Gssat+Cisat+PS100+J100+Vc100+Vo100+Gs100+Ci100+LI_ratio+Narea+Parea+NParea+LMA+Glucose.mg.garea+Fructose.mg.garea+Sucrose.mg.garea+Total.Sugars.mg.garea+Starch.mg.garea+TNC.mg.garea)
summary(pcamodel10)
summary(step(pcamodel10))

pcamodel11 <- lm(Rl~Jsat+LMA+Gs100+Narea)
summary(pcamodel11)
summary(step(pcamodel11))




rda <- rda(respiration,traits)
rda <- rda(traits,respiration)
summary(rda)
plot(rda)

traits <- cbind(avgTleaf,Amax,QY,satPAR,PSsat,Jsat,Vcsat,Vosat,Gssat,Cisat,PS100,J100,Vc100,Vo100,Gs100,Ci100,Rn,LI_ratio,Rlight.Asat,Rdark.Asat,N,P,NP,LMA,Glucose.mg.g,Fructose.mg.g,Sucrose.mg.g,Total.Sugars.mg.g,Starch.mg.g,TNC.mg.g)
traitpca <- princomp(traits)
summary(traitpca,loadings=TRUE,cutoff=.001)
attach(traitpca)
attach(as.data.frame(scores))
plot(traitpca)
biplot(traitpca)
pcascores <- as.data.frame(scores)
attach(pcascores)
plot(traitpca$scores)
write.csv(summary(traitpca,loadings=TRUE,cutoff=.001), "pca.csv")

model <- regsubsets(toolikstand$Rlight.Rdark400~.,data=traits)
summary(model)
plot(model)

model <- gam(Rlight.Rdark400~s(TNC)+s(SLA))
plot(model)










#################################################################################################

#linear regression models



#ambient co2

lm<-lm(Rl.Rn~Amax+QY+PSsat+Jsat+Vcsat+Gssat+Cisat+LMA+Narea+Parea+glucosearea+fructosearea+sucrosearea+starcharea,data=toolik,na.action=na.omit)
summary(lm)
summary(step(lm))
table<-cbind("Rl/Rn",summary(step(lm))$coefficients,glance(step(lm))$adj.r.squared,glance(step(lm))$p.value)


lm2<-lm(Rl~Amax+QY+Vcsat+Cisat+Gssat+LMA+Narea+Parea+glucosearea+fructosearea+sucrosearea+starcharea,data=toolik,na.action=na.omit)
summary(lm2)
summary(step(lm2))
table<-rbind(table,cbind("Rl",summary(step(lm2))$coefficients,glance(step(lm2))$adj.r.squared,glance(step(lm2))$p.value))

lm3<-lm(Rn~Amax+QY+PSsat+Jsat+Vcsat+Gssat+Cisat+LMA+Narea+Parea+glucosearea+fructosearea+sucrosearea+starcharea,data=toolik,na.action=na.omit)
summary(lm3)
summary(step(lm3))
table<-rbind(table,cbind("Rn",summary(step(lm3))$coefficients,glance(step(lm3))$adj.r.squared,glance(step(lm3))$p.value))

lm4<-lm(Rl.Rn~Amaxmass+QY+PSsat+Jsatmass+Vcsatmass+Gssat+Cisat+LMA+N+P+glucose+fructose+sucrose+starch,data=toolik,na.action=na.omit)
summary(lm4)
summary(step(lm4))
table<-rbind(table,cbind("Rl/Rn",summary(step(lm4))$coefficients,glance(step(lm4))$adj.r.squared,glance(step(lm4))$p.value))

lm5<-lm(Rlmass~Amaxmass+QY+PSsat+Jsatmass+Vcsatmass+Gssat+Cisat+LMA+N+P+glucose+fructose+sucrose+starch,data=toolik,na.action=na.omit)
summary(lm5)
summary(step(lm5))
table<-rbind(table,cbind("Rlmass",summary(step(lm5))$coefficients,glance(step(lm5))$adj.r.squared,glance(step(lm5))$p.value))


lm6<-lm(Rnmass~Amaxmass+QY+PSsat+Jsatmass+Vcsatmass+Gssat+Cisat+LMA+N+P+glucose+fructose+sucrose+starch,data=toolik,na.action=na.omit)
summary(lm6)
summary(step(lm6))
table<-rbind(table,cbind("Rnmass",summary(step(lm6))$coefficients,glance(step(lm6))$adj.r.squared,glance(step(lm6))$p.value))
write.csv(table,"lmambient.csv")




#elevated co2


lm<-lm(Rl.Rn~Amax+QY+PSsat+Jsat+Vcsat+Gssat+Cisat+LMA+Narea+Parea+glucosearea+fructosearea+sucrosearea+starcharea,data=elev,na.action=na.omit)
summary(lm)
summary(step(lm))
table<-cbind("Rl/Rn",summary(step(lm))$coefficients,glance(step(lm))$adj.r.squared,glance(step(lm))$p.value)

lm2<-lm(Rl~Amax+QY+Vcsat+Cisat+Gssat+LMA+Narea+Parea+glucosearea+fructosearea+sucrosearea+starcharea,data=elev,na.action=na.omit)
summary(lm2)
summary(step(lm2))
table<-rbind(table,cbind("Rl",summary(step(lm2))$coefficients,glance(step(lm2))$adj.r.squared,glance(step(lm2))$p.value))


lm3<-lm(Rn~Amax+QY+PSsat+Jsat+Vcsat+Gssat+Cisat+LMA+Narea+Parea+glucosearea+fructosearea+sucrosearea+starcharea,data=elev,na.action=na.omit)
summary(lm3)
summary(step(lm3))
table<-rbind(table,cbind("Rn",summary(step(lm3))$coefficients,glance(step(lm3))$adj.r.squared,glance(step(lm3))$p.value))



lm4<-lm(Rl.Rn~Amaxmass+QY+PSsat+Jsatmass+Vcsatmass+Gssat+Cisat+LMA+N+P+glucose+fructose+sucrose+starch,data=elev,na.action=na.omit)
summary(lm4)
summary(step(lm4))
table<-rbind(table,cbind("Rl/Rn",summary(step(lm4))$coefficients,glance(step(lm4))$adj.r.squared,glance(step(lm4))$p.value))


lm5<-lm(Rlmass~Amaxmass+QY+PSsat+Jsatmass+Vcsatmass+Gssat+Cisat+LMA+N+P+glucose+fructose+sucrose+starch,data=elev,na.action=na.omit)
summary(lm5)
summary(step(lm5))
table<-rbind(table,cbind("Rlmass",summary(step(lm5))$coefficients,glance(step(lm5))$adj.r.squared,glance(step(lm5))$p.value))



lm6<-lm(Rnmass~Amaxmass+QY+PSsat+Jsatmass+Vcsatmass+Gssat+Cisat+LMA+N+P+glucose+fructose+sucrose+starch,data=elev,na.action=na.omit)
summary(lm6)
summary(step(lm6))
table<-rbind(table,cbind("Rnmass",summary(step(lm6))$coefficients,glance(step(lm6))$adj.r.squared,glance(step(lm6))$p.value))
write.csv(table,"lmelevated.csv")















lm1<-lm(Rl~Vcsat+LMA+Narea+Parea+glucosearea)
summary(lm1)
summary(step(lm1))





#################################################################################################

#mixed effects model



#ambient co2:

toolik[,9:length(toolik)] <- decostand(toolik[,9:length(toolik)], "standardize")
attach(toolik)

mix<-lme(Rl.Rn~Amax+QY+PSsat+Jsat+Vcsat+Gssat+Cisat+LMA+Narea+Parea+glucosearea+fructosearea+sucrosearea+starcharea,random=~1|species,data=toolik,na.action=na.omit)
summary(mix)
table<-rbind("Rl.Rn",summary(mix)$tTable)


mix2<-lme(Rl~Amax+QY+PSsat+Jsat+Vcsat+Gssat+LMA+Narea+Parea+glucosearea+fructosearea+sucrosearea+starcharea,random=~1|species,data=toolik,na.action=na.omit)
summary(mix2)
table<-rbind(table,"Rl",summary(mix2)$tTable)

mix3<-lme(Rn~Amax+QY+PSsat+Jsat+Vcsat+Gssat+Cisat+LMA+Narea+Parea+glucosearea+fructosearea+sucrosearea+starcharea,random=~1|species,data=toolik,na.action=na.omit)
summary(mix3)
table<-rbind(table,"Rn",summary(mix3)$tTable)




mix4<-lme(Rl.Rn~Amaxmass+QY+PSsat+Jsatmass+Vcsatmass+Gssat+Cisat+LMA+N+P+glucose+fructose+sucrose+starch,random=~1|species,data=toolik,na.action=na.omit)
summary(mix4)
table<-rbind(table,"Rl.Rn",summary(mix4)$tTable)


mix5<-lme(Rlmass~Amaxmass+QY+PSsat+Jsatmass+Vcsatmass+Gssat+Cisat+LMA+N+P+glucose+fructose+sucrose+starch,random=~1|species,data=toolik,na.action=na.omit)
summary(mix5)
table<-rbind(table,"Rl",summary(mix5)$tTable)


mix6<-lme(Rnmass~Amaxmass+QY+PSsat+Jsatmass+Vcsatmass+Gssat+Cisat+LMA+N+P+glucose+fructose+sucrose+starch,random=~1|species,data=toolik,na.action=na.omit)
summary(mix6)
table<-rbind(table,"Rn",summary(mix6)$tTable)
write.csv(table,"mixstand.csv")



#elevated co2:

elev[,9:length(elev)] <- decostand(elev[,9:length(elev)], "standardize")

mix<-lme(Rl.Rn~Amax+QY+PSsat+Jsat+Vcsat+Gssat+Cisat+LMA+Narea+Parea+glucosearea+fructosearea+sucrosearea+starcharea,random=~1|species,data=elev,na.action=na.omit, method="ML")
summary(mix)
table<-rbind("Rl.Rn",summary(mix)$tTable)

mix2<-lme(Rl~Amax+QY+PSsat+Jsat+Vcsat+Gssat+LMA+Narea+Parea+glucosearea+fructosearea+sucrosearea+starcharea,random=~1|species,data= elev,na.action=na.omit)
summary(mix2)
table<-rbind(table,"Rl",summary(mix2)$tTable)

mix3<-lme(Rn~Amax+QY+PSsat+Jsat+Vcsat+Gssat+Cisat+LMA+Narea+Parea+glucosearea+fructosearea+sucrosearea+starcharea,random=~1|species,data= elev,na.action=na.omit)
summary(mix3)
table<-rbind(table,"Rn",summary(mix3)$tTable)




mix4<-lme(Rl.Rn~Amaxmass+QY+PSsat+Jsatmass+Vcsatmass+Gssat+Cisat+LMA+N+P+glucose+fructose+sucrose+starch,random=~1|species,data= elev,na.action=na.omit)
summary(mix4)
table<-rbind(table,"Rl.Rn",summary(mix4)$tTable)


mix5<-lme(Rlmass~Amaxmass+QY+PSsat+Jsatmass+Vcsatmass+Gssat+Cisat+LMA+N+P+glucose+fructose+sucrose+starch,random=~1|species,data= elev,na.action=na.omit)
summary(mix5)
table<-rbind(table,"Rl",summary(mix5)$tTable)


mix6<-lme(Rnmass~Amaxmass+QY+PSsat+Jsatmass+Vcsatmass+Gssat+Cisat+LMA+N+P+glucose+fructose+sucrose+starch,random=~1|species,data= elev,na.action=na.omit)
summary(mix6)
table<-rbind(table,"Rn",summary(mix6)$tTable)
write.csv(table,"mixelevstandardized.csv")









mix<-lmer(Rl.Rn~Amax+QY+PSsat+Jsat+Vcsat+Gssat+Cisat+LMA+Narea+Parea+glucosearea+fructosearea+sucrosearea+starcharea+(1|species),data=toolik)
display(mix)





mixrlrn <- lme(Rl.Rn ~ Amax+QY+PSsat+Jsat+Vcsat+Vosat+Gssat+Cisat+N+P+glucose+fructose+sucrose+starch, data=toolik, random=~1|species, na.action=na.omit)
summary(mix)
anova(mix)
vif(mixrlrn)

mixrl <- lmer(Rl ~ Amax+QY+satPAR+PSsat+Jsat+Vcsat+Vosat+Gssat+Cisat+PS100+J100+Vc100+Vo100+Gs100+Ci100+Rn+LI_ratio+Rlight.Asat+Rdark.Asat+N+P+NP+LMA+Glucose.mg.g+Fructose.mg.g+Sucrose.mg.g+Total.Sugars.mg.g+Starch.mg.g+TNC.mg.g + (1 | avgTleaf), data = toolikstand)
summary(mix)
anova(mix)



#multiple regressions

cors <- abs(cor(toolik[,12:ncol(toolik)]))
rcors <- rcorr(as.matrix(toolik[,12:ncol(toolik)]),type="pearson")
corrgram(toolik[,12:ncol(toolik)])
dev.new(height=12, width=12)
pdf("toolikcorrelations.pdf")
corrgram(cors, upper.panel=panel.pie, lower.panel=NULL, cex.labels=0.1)
dev.off()

model1 <- lm(toolikstand$Rlight.Rdark400~toolikstand$Glucose+toolikstand$Fructose+toolikstand$Sucrose+toolikstand$SLA+toolikstand$Amax400)
stepmodel1 <- step(model1)
plot(model1)
summary(model1)
summary(stepmodel1)
anova(stepmodel1)

model1 <- lm(toolikstand$Rlight400~toolikstand$Glucose+toolikstand$Fructose+toolikstand$Sucrose+toolikstand$SLA+toolikstand$Amax400)
stepmodel1 <- step(model1)
plot(model1)
summary(model1)
summary(stepmodel1)
anova(stepmodel1)

model <- lm(toolikstand$Rlight400~+toolikstand$Fructose+toolikstand$Sucrose+toolikstand$SLA+toolikstand$Amax400)
stepmodel <- step(model)
summary(model)
summary(stepmodel)
anova(stepmodel)

model1 <- lm(toolikstand$Rdark400~toolikstand$Glucose+toolikstand$Fructose+toolikstand$Sucrose+toolikstand$SLA+toolikstand$Amax400)
stepmodel1 <- step(model1)
plot(model1)
summary(model1)
summary(stepmodel1)
anova(stepmodel1)

model2 <- lm(toolik$Rl.Rn~decostand(Amax, "standardize"))
plot(model2)
summary(model2)

model3 <- lm(toolik$Rl~avgTleaf+QY+Rd.actual+Jsat+LMA+N+Amax+P+Glucose.mg.g+Fructose.mg.g+Sucrose.mg.g+Starch.mg.g)
summary(model3)
summary(step(model3))
model3a <- lm(toolik$Rl~LMA)
summary(model3a)

model4 <- lm(toolik$Rn~LMA+N+Amax+P+Glucose.mg.g+Fructose.mg.g+Sucrose.mg.g+Starch.mg.g)
summary(model4)
summary(step(model4))

model5 <- lm(toolik$Rl.Rn~LMA+N+Amax+P+Glucose.mg.g+Fructose.mg.g+Sucrose.mg.g+Starch.mg.g)
summary(model5)
summary(step(model5))

areamodel <- lm(Rl~avgTleaf+Amax+QY+satPAR+PSsat+Jsat+Vcsat+Vosat+Gssat+Cisat+PS100+J100+Vc100+Vo100+Gs100+Ci100+LI_ratio+NAmax+NAsat+PAmax+PAsat+Narea+Parea+NParea+LMA+Glucose.mg.garea+Fructose.mg.garea+Sucrose.mg.garea+Total.Sugars.mg.garea+Starch.mg.garea+TNC.mg.garea)
summary(areamodel)
summary(step(areamodel))
vif(areamodel)

avgTleaf,Amax,QY,satPAR,PSsat,Jsat,Vcsat,Vosat,Gssat,Cisat,PS100,J100,Vc100,Vo100,Gs100,Ci100,LI_ratio,NAmax,NAsat,PAmax,PAsat,Narea,Parea,NParea,LMA,Glucose.mg.garea,Fructose.mg.garea,Sucrose.mg.garea,Total.Sugars.mg.garea,Starch.mg.garea,TNC.mg.garea

areamodelrlrn <- lm(Rl.Rn~avgTleaf+Amax+QY+satPAR+PSsat+Jsat+Vcsat+Vosat+Gssat+Cisat+PS100+J100+Vc100+Vo100+Gs100+Ci100+LI_ratio+Rlight.Asat+Rdark.Asat+#massAmax+massAsat+massJsat+massVcsat+massVosat+massJ100+massVc100+massVo100+massRlight+massRdark+NAmax+NAsat+NRdark+PAmax+PAsat+PRdark+
Narea+Parea+NParea+LMA+Glucose.mg.garea+Fructose.mg.garea+Sucrose.mg.garea+Total.Sugars.mg.garea+Starch.mg.garea+TNC.mg.garea)
summary(areamodel)
summary(step(areamodel))
vif(areamodel)

massmodel <- lm(Rl~avgTleaf+Amax+QY+satPAR+PSsat+Jsat+Vcsat+Vosat+Gssat+Cisat+PS100+J100+Vc100+Vo100+Gs100+Ci100+Rn+LI_ratio+Rlight.Asat+Rdark.Asat+#massAmax+massAsat+massJsat+massVcsat+massVosat+massJ100+massVc100+massVo100+massRlight+massRdark+NAmax+NAsat+NRdark+PAmax+PAsat+PRdark+
Narea+Parea+NParea+LMA+Glucose.mg.garea+Fructose.mg.garea+Sucrose.mg.garea+Total.Sugars.mg.garea+Starch.mg.garea+TNC.mg.garea)

toolik$Rl~avgTleaf+Amax+QY+Rd.actual+satPAR+PSsat+Jsat+Vcsat+Vosat+Gssat+Cisat+PS100+J100+Vc100+Vo100+Gs100+Ci100+Rn+LI_ratio+Rlight.Asat+Rdark.Asat+massAmax+massAsat+massJsat+massVcsat+massVosat+massJ100+massVc100+massVo100+massRlight+massRdark+NAmax+NAsat+NRlight+NRdark+PAmax+PAsat+PRlight+PRdark+N+P+NP+LMA+Glucose.mg.g+Fructose.mg.g+Sucrose.mg.g.+Total.Sugars.mg.g+Starch.mg.g+TNC.mg.g


#multiple regression using mass-based values
model1 <- lm(toolikstand$Rlight.Rdark400~toolikstand$Glucose+toolikstand$Fructose+toolikstand$Sucrose+toolikstand$SLA+toolikstand$amax400w)
stepmodel1 <- step(model1)
plot(model1)
summary(model1)
summary(stepmodel1)
anova(stepmodel1)

model1 <- lm(toolikstand$rl400w~toolikstand$Glucose+toolikstand$Fructose+toolikstand$Sucrose+toolikstand$SLA+toolikstand$amax400w)
stepmodel1 <- step(model1)
plot(model1)
summary(model1)
summary(stepmodel1)
anova(stepmodel1)



plot3d(Rlight.Rdark400, Sucrose, Amax400, col=terrain.colors(length(Rlight.Rdark400)))

model1 <- lm(toolikstand$Rlight.Rdark400~traits$Fructose+traits$Sucrose+traits$Amax400)
stepmodel1 <- step(model1)
plot(model1)
summary(model1)
summary(stepmodel1)
anova(stepmodel1)


model2 <- lm(Rlight.Rdark400~SLA+Amax400)
summary(model2)
plot(model2)

model3 <- lm(toolik400$Rlight.Rdark400~toolik400$Amax400)
summary(model3)
plot(model3)

pcamodel <- lm(toolikstand$Rlight.Rdark400~pcascores$Comp.1*pcascores$Comp.2*pcascores$Comp.3*pcascores$Comp.4*pcascores$Comp.5)
summary(pcamodel)
stepmodel <- step(pcamodel)
summary(stepmodel)
stepmodel$anova

pcamodel1 <- lm(toolikstand$Rlight.Rdark400~pcascores$Comp.1+pcascores$Comp.2+pcascores$Comp.3+pcascores$Comp.4+pcascores$Comp.5)
summary(pcamodel1)
stepmodel1 <- step(pcamodel1)
summary(stepmodel1)
stepmodel1$anova
plot(pcamodel1)

pcamodel2 <- lm(toolikstand$Rlight.Rdark400~traits$Fructose+traits$Glucose+traits$Sucrose+traits$Starch+traits$SLA+traits$Amax400)
summary(pcamodel2)
stepmodel2 <- step(pcamodel2)
summary(stepmodel2)
stepmodel2$anova

traittree <- tree(Rlight.Rdark400~.,data=traits)

solublesugars = fructose + glucose + sucrose

nmds <- metaMDS(traits, dist="euclidean", k=2, trymax=10000, trace=TRUE, noshare=1)
plot(nmds, display=c("sites"))
biplot(nmds)

euclidean <- vegdist(traits, distance = "euclidean")
nmds2<-isoMDS(euclidean,y=cmdscale(euclidean,2),k=2,tol=1e-04,trace=FALSE)
nmds2
plot(nmds2$points[,1],nmds2$points[,2],main="NMDS 2-D Euclidean Distance",xlab=" ",ylab=" ",asp=1,type="b",lwd=2, col=)
stressplot(nmds2,euclidean)
summary(nmds2)

chord <- vegdist(decostand(traits, "norm"), "euclidean")
nmds2<-isoMDS(chord,y=cmdscale(chord,2),k=2,tol=1e-04,trace=FALSE)
nmds2
plot(nmds2$points[,1],nmds2$points[,2],main="NMDS 2-D Chord Distance",xlab=" ",ylab=" ",asp=1,type="b",lwd=2)
stressplot(nmds2,chord)
summary(nmds2)

euclidean <- vegdist(respiration, distance = "euclidean")
nmds2<-isoMDS(euclidean,y=cmdscale(euclidean,2),k=2,tol=1e-04,trace=FALSE)
plot(nmds2$points[,1],nmds2$points[,2],main="NMDS 2-D Euclidean Distance",xlab=" ",ylab=" ",asp=1,type="b",lwd=2)
stressplot(nmds2,euclidean)

#calculate all rates by mass instead of area (i.e. multiply by SLA)
toolik$amax400w <- Amax400*SLA
toolik$rd400w <- Rdark400*SLA
toolik$rl400w <- Rlight400*SLA
toolik$vo400w <- Voforinhib400*SLA
toolik$rd1500w <- Rdark1500*SLA
toolik$rl1500w <- Rlight1500*SLA
toolik$vo1500w <- Voforinhib1500*SLA
toolik$amax1500w <- Amax1500*SLA
attach(toolik)

#redo analyses using mass-based rates
toolik <- aggregate(toolik, list(Species), mean)
toolikw <- cbind(toolik[,5:11], toolik[,13], toolik[,28:31])
names(toolikw)[8] <- "rlrd400"

cor(toolikw)













#partition trait variation
vars <- matrix(NA,2,ncol(alltraits))
colnames(vars) <- names(alltraits)[1:length(names(alltraits))]
rownames(vars) <- c("Species","Within")
for(i in 1:ncol(alltraits)){
	a <- alltraits[,i]
	b <- 100*varcomp(lme(a~1, random=~1|Species, data = alltraits, na.action = na.omit),1)
	vars[1,i] <- b[1]
	vars[2,i] <- b[2]
}
vars <- vars[,order(colnames(vars))]
vars <- vars[order(rownames(vars), decreasing=TRUE),]




par(mar=c(5,5,4,5))
x <- barplot(vars, col=c("darkviolet","deeppink2"), ylab=("proportion of variation (%)"), legend=c("within species","between species"), args.legend=list(x=30,y=110, horiz=TRUE), xaxt="n", space=.5)
text(x, par("usr")[3], labels=paste(colnames(vars)), xpd=TRUE, srt=60, adj=c(1.1,1.1), cex=0.9 )



vpsla <-varcomp(lme(SLA~1, random=~1|Species, data = toolik, na.action = na.omit),1)




#measure phylogenetic signal using Blomberg's K

toolik <- read.csv("nicksummary.csv", header = TRUE)

for (i in 5:ncol(toolik)){
	toolik[,i] <- as.numeric(as.character(toolik[,i]))
}
toolik$Rl.Rn <- toolik$Rl/toolik$Rn
toolik$Rl.Rd400 <- toolik$Rlight400/toolik$Rdark400
toolik$Rl.Rd1500 <- toolik$Rlight1500/toolik$Rdark1500

attach(toolik)


toolik <- aggregate(toolik, list(Species), mean, na.rm=TRUE)
toolik <- toolik[,-c(2:3)]
names(toolik)[1] <- "Species"
toolik$Species <- gsub(" ","_",toolik$Species)
attach(toolik)

#write.csv(toolik$Group.1,"phylomatic.csv")

tree <- read.tree("tree2.txt")
setwd("/Users/Nick/Dropbox/Work/Columbia/Spring 2008/Alaska/phylogeny")
tree <- read.tree("molecular_tree.new")
setwd("/Users/Nick/Dropbox/Work/Columbia/Spring 2008/Alaska/phylogeny/phylocom/run1")
tree <- read.tree("prior.rtf")
tree <- scan("prior.rtf",what="string")
treezy <- substr(tree[10], 1, nchar(tree[10])-1)


par(bg = 'black')
plot(tree, edge.color='white', tip.color='white', cex=0.9)

tootrix <- data.matrix(toolik)
row.names(tootrix) <- Species
tooliktree <- match.phylo.comm(tree,t(tootrix[,-1]))

#toolikcomm <- dcast(toolik,Species ~ variable, mean)

toolik[,1] <- sort(tree$tip.label)
toolikord <- toolik[match(tree$tip.label, toolik[,1]),]
attach(toolikord)

#levels(toolik[,1])[39] = "vaccinium_vitis-idae"


Rlightstand <- Rlight400 + 10
rlrdstand <- Rlight.Rdark400 + 10
rl400wstand <- rl400w + 800
rl <- Rl + 10

Krand(tooliktree$phy,TNC.mg.g,n=1000)
Krand(tooliktree$phy,Soluble,n=1000)
Krand(tooliktree$phy,Fructose.mg.g,n=1000)
Krand(tooliktree$phy,Glucose,n=1000)
Krand(tooliktree$phy,Sucrose,n=1000)
Krand(tooliktree$phy,Starch,n=1000)
Krand(tooliktree$phy,LMA,n=1000)
Krand(tooliktree$phy,Amax,n=10)
k <- Krand(tooliktree$phy,toolikord$Rl.Rd400,n=1000)
k <- Krand(tooliktree$phy,toolikord$Rl.Rd1500,n=1000)

k <- Krand(tooliktree$phy,Rl.Rn,n=10)
Krand(tree,Rdark400,n=1000)
Krand(tree,Rlightstand,n=1000)
Krand(tree,Voforinhib400,n=1000)
Krand(tree,amax400w,n=1000)
Krand(tree,rd400w,n=1000)
Krand(tree,rl400w,n=1000)
Krand(tree,vo400w,n=1000)
Krand(tree,rlrdstand,n=1000)
Krand(tree,rl400wstand,n=1000)
Krand(tree,N,n=1000)
rl <- as.numeric(decostand(Rl,"range"))
rl2 <- as.numeric(scale(Rl))
Krand(tree,log(rl),n=10)
Krand(tree,abs(Rn),n=10)


nam <- function(x){
	for(i in 1:length(x)){
	if(x[i]== <- mean(x, na.rm=TRUE)		
	}
}

# ks <- matrix(NA,55,6)
# for(i in 10:ncol(toolikord)){
	# ks[i-9,] <- c(names(toolikord[i]),as.numeric(Krand(tree,toolikord[,i],n=1000)))
# }

ks <- matrix(NA,65,7)
for(i in 9:nrow(tooliktree$comm)){
	ks[i-8,] <- c(rownames(tooliktree$comm)[i],as.numeric(Krand(tooliktree$phy,as.numeric(tooliktree$comm[i,]),n=10000)),NA)
	ks[i-8,7] <- ifelse(ks[i-8,4] < ks[i-8,2] & ks[i-8,2] < ks[i-8,5],"","*")
}
colnames(ks) <- c("trait","K.obs","K.null.mean","K.null.lwr","K.null.upr","n.rand","significant")
ks <- as.data.frame(ks)
ks$K.obs <- as.numeric(as.character(ks$K.obs))
ks$K.null.lwr <- as.numeric(as.character(ks$K.null.lwr))
ks$K.null.upr <- as.numeric(as.character(ks$K.null.upr))
attach(ks)
ks <- as.data.frame(ks[order(ks$trait),])
ggplot(ks, aes(x=K.obs, y=trait)) + 
geom_point(aes(color=ifelse(ks$significant=="*","blue","black"))) +
geom_errorbarh(aes(xmin=K.null.lwr, xmax=K.null.upr),width=0.25, color=ifelse(ks$significant=="*","blue","black")) +
theme(axis.text=element_text(size=7)) +
theme(axis.text.y=element_text(face=ifelse(ks$significant=="*","bold","plain"), color=ifelse(ks$significant=="*","blue","black"))) + 
scale_x_log10() + xlab("log Blomberg's K") +
scale_color_manual(values=c("black","blue"),labels=c("p ≥ 0.05","p < 0.05")) +
theme(legend.title=element_blank())
ggsave("ks plot.tiff", path="/Users/Nick/Dropbox/Work/Toronto/Thesis/figures")
setwd("/Users/Nick/Dropbox/Work/Toronto/Thesis/figures")
write.csv(ks,"ks table.csv")

Kcalc(tree,matrix(log(decostand(toolikord[,i],"range")+0.00000001

#phylogenetically independent contrasts

pictree <- multi2di(tooliktree$phy)
ps <- matrix(NA,65,6)
for(i in 9:nrow(tooliktree$comm)){
	ps[i-8,] <- c(rownames(tooliktree$comm)[i],as.numeric(pfunc(as.numeric(tooliktree$comm[i,]), pictree,n=10000)),NA)
	ps[i-8,6] <- ifelse(ps[i-8,4] < ps[i-8,2] & ps[i-8,2] < ps[i-8,5],"","*")
}
colnames(ps) <- c("trait","pic.obs","pic.null.mean","pic.null.05","pic.null.95","significant")
ps <- as.data.frame(ps)
ps$pic.obs <- as.numeric(as.character(ps$pic.obs))
ps$pic.null.05 <- as.numeric(as.character(ps$pic.null.05))
ps$pic.null.95 <- as.numeric(as.character(ps$pic.null.95))
attach(ps)
ps <- as.data.frame(ps[order(ps$trait),])
ggplot(ps, aes(x=pic.obs, y=trait)) + 
geom_point(aes(color=ifelse(ps$significant=="*","blue","black"))) +
geom_errorbarh(aes(xmin=pic.null.05, xmax=pic.null.95),width=0.25, color=ifelse(ps$significant=="*","blue","black")) +
theme(axis.text=element_text(size=7)) +
theme(axis.text.y=element_text(face=ifelse(ps$significant=="*","bold","plain"), color=ifelse(ps$significant=="*","blue","black"))) + 
scale_x_log10() + xlab("log phylogenetically independent contrasts") +
scale_color_manual(values=c("black","blue"),labels=c("p ≥ 0.05","p < 0.05")) +
theme(legend.title=element_blank())
ggsave("ps plot.tiff", path="/Users/Nick/Dropbox/Work/Toronto/Thesis/figures")
setwd("/Users/Nick/Dropbox/Work/Toronto/Thesis/figures")
write.csv(ps, "ps table.csv")

#bivariate regressions

plot(N, Rl)
abline(lm(Rl ~ N))
summary(lm(Rl ~ N))

plot(P, Rl)
abline(lm(Rl ~ P))
summary(lm(Rl ~ P))

plot(Amax, Rl)
abline(lm(Rl ~ Amax))
summary(lm(Rl ~ Amax))

plot(Jsat, Rl)
abline(lm(Rl ~ Jsat))
summary(lm(Rl ~ Jsat))

plot(Vcsat, Rl)
abline(lm(Rl ~ Vcsat))
summary(lm(Rl ~ Vcsat))

plot(Vcsat, Jsat)
abline(lm(Jsat ~ Vcsat))
summary(lm(Jsat ~ Vcsat))

plot(Starch.mg.garea, Rl)
abline(lm(Rl ~ Starch.mg.garea))
summary(lm(Rl ~ Starch.mg.garea))

plot(Starch.mg.garea, Rn)
abline(lm(Rn ~ Starch.mg.garea))
summary(lm(Rn ~ Starch.mg.garea))

plot(N, Rl.Rn)
abline(lm(Rl.Rn ~ N))
summary(lm(Rl.Rn ~ N))


plot(amax400w, rl400w)
abline(lm(rl400w ~ amax400w))
text(3500,-100,paste("r2 = ", round(summary(lm(rl400w ~ amax400w))$adj.r.squared, digits=2)))
summary(lm(rl400w ~ amax400w))

plot(SLA,Rlight400, xlab=expression(paste(SLA," ",(m^2/g))), ylab=expression(paste(R[L]," [",mu,mol,"/(",m^2,"s)]")))
abline(lm(Rlight400 ~ SLA))
rsq <- round(summary(lm(Rlight400 ~ SLA))$adj.r.squared, digits=2)
text(250,-8,expression(r^2==rsq))
summary(lm(Rlight400 ~ SLA))

plot(SLA,Rdark400)
abline(lm(Rdark400 ~ SLA))
summary(lm(Rdark400 ~ SLA))

plot(SLA, Rlight.Rdark1500)
abline(lm(Rlight.Rdark1500 ~ SLA))
summary(lm(Rlight.Rdark1500 ~ SLA))

plot(SLA,Rlight1500)
abline(lm(Rlight1500 ~ SLA))
summary(lm(Rlight1500 ~ SLA))

plot(SLA,Rdark1500)
abline(lm(Rdark1500 ~ SLA))
summary(lm(Rdark1500 ~ SLA))

plot(Rlight.Rdark400 ~ Amax400)
abline(lm(Rlight.Rdark400 ~ Amax400))
summary(lm(Rlight.Rdark400 ~ Amax400))

plot(Voforinhib400, Rlight.Rdark400)
abline(lm(Rlight.Rdark400 ~ Voforinhib400))
summary(lm(Rlight.Rdark400 ~ Voforinhib400))

Voforinhib400 <- rm.outlier(Voforinhib400, fill = TRUE)

plot(Voforinhib1500, Rlight.Rdark1500)
abline(lm(Rlight.Rdark1500 ~ Voforinhib1500))
summary(lm(Rlight.Rdark1500 ~ Voforinhib1500))

plot(Sucrose.mg.g,Rl)
abline(lm(Rl~Sucrose.mg.g))
summary(lm(Rl~Sucrose.mg.g))

plot(Sucrose.mg.g,Rl.Rn)
abline(lm(Rl.Rn~Sucrose.mg.g))
summary(lm(Rl.Rn~Sucrose.mg.g))

plot(LMA,Rl)
abline(lm(Rl~LMA))
summary(lm(Rl~LMA))

#2009

attach(toolik)
toolik <- aggregate(toolik, list(Species), mean)
toolik$rnew <- toolik$Rlight.Rdark400/toolik$Amax400
attach(toolik)
head(toolik)

mahmoud <- read.csv("mahmoud2500XT.csv")
attach(mahmoud)
head(mahmoud)

#create species list for Imnavait
write.csv(unique(mahmoud$Species),"imnavaitspecies.csv")

cover <- cbind(data.frame(mahmoud[,3],mahmoud[,6],mahmoud[,10]))
names(cover)[1] <- "transect"
names(cover)[2] <- "species"
names(cover)[3] <- "cover"
attach(cover)
head(cover)

transectmeans <- aggregate(cover, by=list(Transect,Species), mean)
transectmeans <- cbind(data.frame(transectmeans[,1],transectmeans[,2],transectmeans[,5]))
names(transectmeans)[1] <- "transect"
names(transectmeans)[2] <- "species"
names(transectmeans)[3] <- "cover"
transectmeans <- as.matrix(transectmeans)

barplot(t(transectmeans), beside=FALSE)

write.csv(transectmeans, "transectmeans.csv")

toolik <- read.csv("tooliknumbers.csv", header = TRUE)
attach(toolik)
head(toolik)
toolik <- aggregate(toolik, list(Species), mean)
rdark <- cbind(data.frame(toolik[,1],toolik[,14]))
names(rdark)[1] <- "species"
names(rdark)[2] <- "rdark"

darkcover <- merge(transectmeans, rdark, by.x="species", by.y="species")
darkcover$rweighted <- darkcover$cover*darkcover$rdark
darktransects <- cbind(data.frame(darkcover[,2],darkcover[,5]))
names(darktransects)[1] <- "transect"
names(darktransects)[2] <- "rweighted"
attach(darktransects)

darktransects <- tapply(darktransects$rweighted,darktransects$transect,sum)
write.csv(darktransects, "transectrates.csv")

#calculate transect rlight
rlight <- cbind(data.frame(toolik[,1],toolik[,15]))
names(rlight)[1] <- "species"
names(rlight)[2] <- "rlight"
lightcover <- merge(transectmeans, rlight, by.x="species", by.y="species")
lightcover$rweighted <- lightcover$cover*lightcover$rlight
lighttransects <- cbind(data.frame(lightcover[,2],lightcover[,5]))
names(lighttransects)[1] <- "transect"
names(lighttransects)[2] <- "rweighted"
attach(lighttransects)

lighttransects <- tapply(lighttransects$rweighted,lighttransects$transect,sum)
write.csv(lighttransects, "lighttransectrates.csv")



xyplot(Area..LAM. ~ Area..PC. | Code)
summary(xyplot(Area..LAM. ~ Area..PC. | Code))

areas <- as.data.frame(cbind(Area.LAM[1:190],Area.PC[1:190]))
areas <- rm.outlier(areas)
attach(areas)

plot(Area.PC[1:190], Area.LAM[1:190],xlab=expression("Area from percent cover "~(cm^2)),ylab=expression("Area from leaf area meter "~(cm^2)))
abline(lm(Area.LAM[1:190] ~ Area.PC[1:190]))
summary(lm(areas$V1 ~ areas$V2))
plot(areas$V2,areas$V1)

plot(Percent.Cover.Top, Percent.Cover.LAM)
abline(lm(Percent.Cover.LAM ~ Percent.Cover.Top))
summary(lm(Percent.Cover.LAM ~ Percent.Cover.Top))

leaves <- read.csv("D98 Small envelopes 1.csv")
unique(leaves[,2])

Krand<- function (phy,com,n.rand=100,sub=FALSE,fig=TRUE) {
	
	if (sub==TRUE) {
		source("/Volumes/cadotte/Documents/R/phylogeny code/subTree.R")
		phy<-subTree(phy,com[,3])
		}
		
	trait <- com
	
	#variable has to be a matrix
	trait <- matrix(trait,ncol=1)
		
	K.obs<-Kcalc(trait,phy)
	K.obs<-data.frame(K.obs)	
	
	kvalues<-numeric(length(n.rand))
	
	for (i in 1:n.rand) {
 		kvalues[i]<-Kcalc(sample(trait),phy)
		}

	#null stats
	K.null<-kvalues
	K.null.mean<-mean(K.null)
	K.null.95<-quantile(K.null,probs=c(0.025,0.975))
	names(K.null.95)<-NULL
	
	results<-cbind(K.obs,K.null.mean,
		K.null.lwr=K.null.95[1],
		K.null.upr=K.null.95[2],
		n.rand)
					
	if (fig==TRUE) {
		par(cex.lab=1.5, cex.axis=1.2, cex.main=1.5, oma=c(0,1.5,0,0))
		hist(K.null, xlab="K values",
		main="Null K values; dashed line = K obs; dotted = 95% CI")
		abline(v=K.obs, lty="dashed")
		abline(v=K.null.95[1],lty="dotted")
		abline(v=K.null.95[2],lty="dotted")
		}	
	
	return(results)	
	}


####function to calculate the sum of absolute phylogentically independent contrasts values and compare them to null values. Requires the package "ape". Arguments: dat is a vector of values sorted according to species in phylogeny, tree is the phylogeny and n.rand is the number of randomizations. Example sorting: dat<-trait.data$trait1[match(tree$tip.label,trait.data$species)]

pfunc<-function(dat,tree,n.rand=1000){
	
		pic.null<-NULL

	for (i in 1:n.rand){
		tmp<-sample(dat,size=length(dat))
		pic.null[i]<-sum(abs(pic(tmp,tree)))
		}

	results<-data.frame(pic.obs=sum(abs(pic(dat,tree))),
		pic.null.mean=mean(pic.null),
		pic.null.05=quantile(pic.null,probs=c(0.05,0.95))[1],
		pic.null.95=quantile(pic.null,probs=c(0.05,0.95))[2],
		row.names=NULL
		)
	
	return(results)
	}




	
#partition variance in SLA and Amax

nitrogen <- read.csv("Nitrogen.csv", header = TRUE)
attach(nitrogen)

varcomp.SLA<-varcomp(lme(LMA~1, random=~1|Species, data = toolik, na.action = na.omit),1)
varcomp.SLA
variation <- rbind(varcomp.SLA)

varcomp.amax<-varcomp(lme(Amax~1, random=~1|Species, data = toolik, na.action = na.omit),1)
varcomp.amax
variation <- rbind(variation,varcomp.amax)

varcomp.Rl<-varcomp(lme(Rl~1, random=~1|Species, data = toolik, na.action = na.omit),1)
varcomp.Rl
variation <- rbind(variation,varcomp.Rl)

varcomp.Rn<-varcomp(lme(Rn~1, random=~1|Species, data = toolik, na.action = na.omit),1)
varcomp.Rn
variation <- rbind(variation,varcomp.Rn)

varcomp.Rl.Rn<-varcomp(lme(Rl.Rn~1, random=~1|Species, data = toolik, na.action = na.omit),1)
varcomp.Rl.Rn
variation <- rbind(variation,varcomp.Rl.Rn)

varcomp.N<-varcomp(lme(N~1, random=~1|Species, data = nitrogen[(1&41:58),], na.action = na.omit),1)
varcomp.N
variation <- rbind(variation,varcomp.N)

varcomp.P<-varcomp(lme(P~1, random=~1|Species, data = nitrogen[(1&41:58),], na.action = na.omit),1)
varcomp.P
variation <- rbind(variation,varcomp.P)

varcomp.massAmax<-varcomp(lme(massAmax~1, random=~1|Species, data = toolik, na.action = na.omit),1)
varcomp.massAmax
variation <- rbind(variation,varcomp.massAmax)

varcomp.massRlight<-varcomp(lme(massRlight~1, random=~1|Species, data = toolik, na.action = na.omit),1)
varcomp.massRlight
variation <- rbind(variation,varcomp.massRlight)

varcomp.massRdark<-varcomp(lme(massRdark~1, random=~1|Species, data = toolik, na.action = na.omit),1)
varcomp.massRdark
variation <- rbind(variation,varcomp.massRdark)

varcomp.NAmax<-varcomp(lme(NAmax~1, random=~1|Species, data = toolik, na.action = na.omit),1)
varcomp.NAmax
variation <- rbind(variation,varcomp.NAmax)

varcomp.NRlight<-varcomp(lme(NRlight~1, random=~1|Species, data = toolik, na.action = na.omit),1)
varcomp.NRlight
variation <- rbind(variation,varcomp.NRlight)

varcomp.NRdark<-varcomp(lme(NRdark~1, random=~1|Species, data = toolik, na.action = na.omit),1)
varcomp.NRdark
variation <- rbind(variation,varcomp.NRdark)

varcomp.PAmax<-varcomp(lme(PAmax~1, random=~1|Species, data = toolik, na.action = na.omit),1)
varcomp.PAmax
variation <- rbind(variation,varcomp.PAmax)

varcomp.PRlight<-varcomp(lme(PRlight~1, random=~1|Species, data = toolik, na.action = na.omit),1)
varcomp.PRlight
variation <- rbind(variation,varcomp.PRlight)

varcomp.PRdark<-varcomp(lme(PRdark~1, random=~1|Species, data = toolik, na.action = na.omit),1)
varcomp.PRdark
variation <- rbind(variation,varcomp.PRdark)

varcomp.NP<-varcomp(lme(NP~1, random=~1|Species, data = toolik, na.action = na.omit),1)
varcomp.NP
variation <- rbind(variation,varcomp.NP)

variation <- as.data.frame(variation)
attach(variation)
write.csv(variation,"variancepartitioning.csv")

sugars <- as.data.frame(cbind(Species,Glucose.mg.g,Fructose.mg.g,Sucrose.mg.g,Total.Sugars.mg.g,Starch.mg.g))
attach(sugars)
sugarsordered <- as.data.frame(sugars[order(sugars$Species),])


toolikft <- cbind(data.frame(toolik[,1], toolik[,10:11]))
names(toolikft)[1] <- "species"
head(toolikft)
toolikfd <- dist(toolikft)
head(toolikfd)

counts <- cbind(data.frame(toolik[,1], toolik[,6]))
abun <- aggregate(counts, list(Species), FUN = sum)
abun <- cbind(matrix(abun[,1], abun[,3]))
names(abun)[1] <- "species"
tooliksp <- cbind(data.frame(toolik[,1], toolik[,11:12]))
names(tooliksp)[1] <- "species"
head(tooliksp)
toolikd <- dist(tooliksp)

FDin(toolikfd)
FDis(toolikfd, abun)

