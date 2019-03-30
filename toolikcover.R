#this takes the 2006lgcover file and averages the cover of each species from the various control plots


setwd("/Users/Nick/Documents/School Work/Columbia/Spring 2008/Alaska")

ak <- read.csv("2006lgcover.csv")
names(ak)
attach(ak)

CT = subset(ak, Treatment == "CT")

rel <- tapply(CT$Relative, CT$Species, mean)

splist <- sort(rel, decreasing = TRUE)

write.csv(splist, "splist.csv")