dat<-read.csv("AvC14_geno_trans_pre_pos.csv",header=TRUE)
#if 8 variables
dat<-dat[,-2]

uniqlist<-unique(dat[,1])
#local
done<-matrix(0,nrow=1,ncol=7)
colnames(done)<-colnames(dat)
for (i in 1:length(uniqlist)){
  temp<-dat[dat[,1]==uniqlist[i],]
  chrom<-as.character(temp[1,4])
  pos.min<-temp[1,6]-25001
  pos.max<-temp[1,7]+25001
  temp1<-temp[temp[,2]==chrom & temp[,3]>=pos.min & temp[,3]<=pos.max,]
  done<-rbind(done,temp1)
}
#long range
done<-matrix(0,nrow=1,ncol=7)
colnames(done)<-colnames(dat)
for (i in 1:length(uniqlist)){
  temp<-dat[dat[,1]==uniqlist[i],]
  chrom<-as.character(temp[1,4])
  pos.min<-temp[1,6]-150001
  pos.max<-temp[1,7]+150001
  temp1<-temp[temp[,2]==chrom & temp[,3]>=pos.min & temp[,3]<=pos.max,]
  done<-rbind(done,temp1)
}

