setwd("~/Documents/R Stuff/Transcriptomics FLAM/75trim")
exp<-read.csv(file="AvC14_genelist.csv",header=TRUE)
exp.de<-read.csv(file="fecun_d14_alltime.csv",header=TRUE)
exp.de1<-exp[(exp[,1] %in% exp.de[,2]),]
exp.de2<-exp.de1[,2:21]
allele.freqs<-t(exp.de2)
loci.list<-as.character(exp.de1[,1])
colnames(allele.freqs)<-loci.list
bcd.50kb.sparse.matrix.af<-allele.freqs
#sim.pheno.all<-read.csv("mortality_sorted.csv",header=TRUE)
sim.pheno.all<-read.csv("fecundity_sorted.csv",header=TRUE)
a.list<-(1:19)*5/100
#looped
corr.list<-NULL
conf.list<-matrix("q",nrow=8,ncol=2)
dumb<-sapply(1:(length(sim.pheno.all[1,])-1),function(y){
  #sim.pheno<-log(sim.pheno.all[,y+1])
  sim.pheno<-sim.pheno.all[,y+1]
  library(flam)
  bcd50.GGD.error.af<- sapply(1:length(a.list),function(x){
    ggd.flam<- flamCV(bcd.50kb.sparse.matrix.af,sim.pheno,alpha=a.list[x],n.fold=5,seed=1,method="GGD")
    best.hat<- ggd.flam$index.cv
    ggd.flam$mean.cv.error[best.hat]
  })
  a.hat.af<- which.min(bcd50.GGD.error.af)
  bcd50.GGD.flam.af<- flamCV(bcd.50kb.sparse.matrix.af,sim.pheno,alpha=a.list[a.hat.af],n.fold=5,seed=1,method="GGD")
  a.pred<- NULL # predicted A phenotypes
  c.pred<- NULL # predicted C phenotypes
  for (i in 1:5){
    temp<- c((i-1)*2+1,(i-1)*2+2,(i-1)*2+11,(i-1)*2+12) # This generates the four populations to exclude, 2A's and 2C's
    train.mtx<- bcd.50kb.sparse.matrix.af[-temp,]
    train.pheno<- sim.pheno[-temp]
    ggd.flam<- flam(train.mtx,train.pheno,alpha=bcd50.GGD.flam.af$alpha,method="GGD")
    a.c.pred<- predict(ggd.flam,new.x=bcd.50kb.sparse.matrix.af[temp,],lambda=bcd50.GGD.flam.af$lambda.cv,alpha=bcd50.GGD.flam.af$alpha)
    a.pred<- c(a.pred,a.c.pred[1:2])
    c.pred<- c(c.pred,a.c.pred[3:4])
  }
  corr.list<<-c(corr.list,cor(c(a.pred,c.pred),sim.pheno))
  conf<-cor.test(c(a.pred,c.pred),sim.pheno)
  conf.list[y,1]<<-conf$conf.int[1]
  conf.list[y,2]<<-conf$conf.int[2]
})
corr.fin<-cbind(colnames(sim.pheno.all)[2:9],corr.list)

sim.pheno<-log(sim.pheno.all[,4])
#flam
library(flam)
bcd50.GGD.error.af<- sapply(1:length(a.list),function(x){
  ggd.flam<- flamCV(bcd.50kb.sparse.matrix.af,sim.pheno,alpha=a.list[x],n.fold=5,seed=1,method="GGD")
  best.hat<- ggd.flam$index.cv
  ggd.flam$mean.cv.error[best.hat]
})
a.hat.af<- which.min(bcd50.GGD.error.af)
bcd50.GGD.flam.af<- flamCV(bcd.50kb.sparse.matrix.af,sim.pheno,alpha=a.list[a.hat.af],n.fold=5,seed=1,method="GGD")
a.pred<- NULL # predicted A phenotypes
c.pred<- NULL # predicted C phenotypes
for (i in 1:5){
  temp<- c((i-1)*2+1,(i-1)*2+2,(i-1)*2+11,(i-1)*2+12) # This generates the four populations to exclude, 2A's and 2C's
  train.mtx<- bcd.50kb.sparse.matrix.af[-temp,]
  train.pheno<- sim.pheno[-temp]
  ggd.flam<- flam(train.mtx,train.pheno,alpha=bcd50.GGD.flam.af$alpha,method="GGD")
  a.c.pred<- predict(ggd.flam,new.x=bcd.50kb.sparse.matrix.af[temp,],lambda=bcd50.GGD.flam.af$lambda.cv,alpha=bcd50.GGD.flam.af$alpha)
  a.pred<- c(a.pred,a.c.pred[1:2])
  c.pred<- c(c.pred,a.c.pred[3:4])
}
corr.list<-NULL
for (x in 2:9){
  corr.list<-c(corr.list,cor(c(a.pred,c.pred),log(sim.pheno.all[,x])))
}
#for (x in 2:9){
  #corr.list<-c(corr.list,cor(c(a.pred,c.pred),sim.pheno.all[,x]))
#}
corr.fin<-cbind(colnames(sim.pheno.all)[2:9],corr.list)
