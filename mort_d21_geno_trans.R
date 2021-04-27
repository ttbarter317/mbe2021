setwd("~/Dropbox/Research Project/FLAM/Geno vs Trans")
trans<-read.csv("mort_d21_alltime.csv",header=TRUE)
trans<-as.matrix(trans[1:19,3])
snp<-read.csv("mort_d21_geno_only.csv",header=TRUE)
together<-read.csv("geno_trans_snp_causal.csv",header=TRUE)
trans.together<-together[together[,1] %in% trans,]
library(tidyr)
