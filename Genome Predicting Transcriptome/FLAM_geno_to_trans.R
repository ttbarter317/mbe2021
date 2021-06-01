setwd("~/Documents/R Stuff/Transcriptomics FLAM/75trim")
# FLAM for EXPERIMENTAL EVOLUTION
# The genetic data consists of a matrix with n-rows and m-columns
# The allele frequencies must be inserted into the matrix "allele.freqs"
# There also needs to be a variable "loci.list" with the names of each locus
# in the same order as the appear in "allele.freqs" so after permuting
# columns we can reconstruct the original loci names.

# In the code below coverage values are read in and converted to allele frequencies.
snp.data<- read.table(file ="A_C_sigSNPs_50kb_v3_pvals.txt",header=T)
# Make a list of the loci designations
loci.list<- sapply(1:length(snp.data[,1]),function(x){
  return(paste(snp.data[x,1],"_",snp.data[x,2]))
})
# Now calculate a matrix of allele frequencies. Each row is a single locus with 20 frequencies
# the dim(allele.freqs) is npops by nloci
col.num<- length(snp.data[1,])-1
allele.freqs<- sapply(1:length(snp.data[,1]),function(x){
  temp<- snp.data[x,5:col.num]
  temp1<- sapply(1:20, function(i) temp[2*i-1]/(temp[2*i-1]+temp[2*i]))
  return(as.numeric(temp1))
})
colnames(allele.freqs)<- loci.list

######################################
# Next read in the phenotype data.   #
# In this case I use only one column #
# of a data frame with 4 columns     #
######################################

#phenotype
exp<-read.csv(file="AvC21_genelist.csv",header=TRUE)
exp.de<-read.csv(file="AvC21_75trim_lme_de5_genes.csv",header=TRUE)
exp.de1<-exp[(exp[,1] %in% exp.de[,1]),]
exp.de2<-exp.de1[,2:21]
pheno<-t(exp.de2)
loci.list1<-as.character(exp.de1[,1])
colnames(pheno)<-loci.list1
#pheno.data<-log(pheno.data)

# Now we are ready to start the FLAM analysis
library(foreach)
library(doParallel)
registerDoParallel(cores=40) # This should be set to the machine capabiliities 
flam.exp.evo.version<- version$version.string
Nper<- 100 # This is the number of times we will permute the columns of the genetic matrix
# If possible I would consider jackking this up to perhpas 500
set.seed(100)
a.list<-(1:19)*5/100 # This is the grid of a values and we will search for the minimum CV error 
# over all 19 values of a
# If there is an error in any of the loops of foreach the error message
# will be passed on to per.list but the other loops will continue going
fin<-matrix(0,nrow=length(pheno[1,]),ncol=3)
fin[,1]<-loci.list1
dumb<-sapply(1:length(pheno[1,]),function(q){
  pheno.data<-pheno[,q]
  per.list<- foreach(k=1:Nper,.errorhandling = c('pass')) %dopar%{ 
    library(flam)
    loci.per<- sample(1:length(allele.freqs[1,]),length(allele.freqs[1,])) #permute the columns of the genetic matrix
    sim.gen<- allele.freqs[,loci.per]
    a.grid<- list() 
    for (x in 1:length(a.list)){
      sim.flamCV<- flamCV(sim.gen,pheno.data,alpha=a.list[x],n.fold=5,seed=1,method="BCD")#we need to use the same seed for each a
      #Output list with sparse loci and CV error
      best.hat<- sim.flamCV$index.cv
      best.sparse<- sim.flamCV$flam.out$non.sparse.list[[best.hat]]
      best.a<- sim.flamCV$mean.cv.error[best.hat]
      a.grid[[x]]<- list(best.a,loci.list[loci.per[best.sparse]])
    }
    a.min<- NULL
    for (i in 1:length(a.grid)) a.min<- c(a.min,unlist(a.grid[[i]][1]))
    best.grid<- which.min(a.min)
    a.grid[[best.grid]][2]
  }
  # Now search through the 100 simulations and construct a list of all loci that appear on any list
  loci.counts<- NULL
  for (i in 1:Nper) loci.counts<- c(loci.counts,unlist(per.list[[i]]))
  u.loci<- sort(unique(loci.counts))
  # Now get a frequency count of the unique loci
  u.loci.freq<- NULL
  for (i in 1:length(u.loci)) u.loci.freq<- c(u.loci.freq,length(loci.counts[loci.counts==u.loci[i]]))
  sparse.list.counts<- cbind(u.loci,u.loci.freq)
  if (length(sparse.list.counts[1,])>=1.99){
    crit.k<- 0.5 # this will create a list of all SNPs satisfying a 50% criteria
    max.freq<- max(as.numeric(sparse.list.counts[,2]))*crit.k
    reduced.sparse.list<- sparse.list.counts[as.numeric(sparse.list.counts[,2])>max.freq,1]
    fin[q,2]<<-length(as.matrix(reduced.sparse.list)[,1])
    fin[q,3]<<-paste0(as.matrix(reduced.sparse.list)[,1],collapse=";")
  }
})
