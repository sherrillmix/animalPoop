library(dnar)
library(parallel)
fastqs<-list.files('split','fastq.gz',full.names=TRUE)
readCounts<-unlist(mclapply(fastqs,function(xx){
  nrow(read.fastq(xx))
},mc.cores=30))
names(readCounts)<-fastqs
pairs<-sub('_[12].fastq.gz$','',basename(fastqs))
if(any(table(pairs)!=2))stop("Messed up left right pair")
if(any(tapply(readCounts,pairs,function(x)all(x!=x[1]))))stop('Count disagreement')
pairCounts<-tapply(readCounts,pairs,'[[',1)

write.csv(data.frame('sample'=pairs,'file'=basename(fastqs),'count'=readCounts),'counts.csv',row.names=FALSE)

samples<-read.csv('Scott_Island_Biogeo_Combined_MFs_V5.csv')
samples$readCounts<-pairCounts[samples$X.SampleID]
samples$qpcr<-apply(sample[,c('X16S.qPCR.copies.per.reaction.replicate.1','X16S.qPCR.copies.per.reaction.replicate.2')],1,mean)

table(samples$readCounts>100,samples$PlateNumber)

