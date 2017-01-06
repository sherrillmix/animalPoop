library(dnar)
library(parallel)
fastqs<-c(list.files('split','fastq.gz',full.names=TRUE),list.files('split2','fastq.gz',full.names=TRUE))
readCounts<-unlist(mclapply(fastqs,function(xx){
  nrow(read.fastq(xx))
},mc.cores=30))
names(readCounts)<-fastqs
pairs<-sub('_[12].fastq.gz$','',basename(fastqs))
runs<-ifelse(dirname(fastqs)=='split',1,ifelse(dirname(fastqs)=='split2',2,NA))
if(any(!table(pairs,runs)%in% c(0,2)))stop("Messed up left right pair")
if(any(tapply(readCounts,paste(runs,pairs),function(x)all(x!=x[1]))))stop('Count disagreement')
pairCounts<-tapply(readCounts,list(pairs,sprintf('run%d',runs)),'[[',1)

write.csv(pairCounts,'counts.csv')

samples<-read.csv('sherrill-MixLauder_islandGut.csv',stringsAsFactors=FALSE)
samples$readCounts<-apply(pairCounts,1,max,na.rm=TRUE)[samples$X.SampleID]
samples$qpcr<-apply(samples[,c('X16S.qPCR.copies.per.reaction.replicate.1','X16S.qPCR.copies.per.reaction.replicate.2')],1,mean)

table(samples$readCounts>100,samples$PlateNumber)

i1<-read.fastq('data/Undetermined_S0_L001_I1_001.fastq.gz')
i2<-read.fastq('data/Undetermined_S0_L001_I2_001.fastq.gz')
bars<-paste(i1$seq,i2$seq)
tail(sort(table(bars[!bars %in% paste(samples$Index1Sequence,samples$Index2Sequence)])))

tab1<-sort(table(i1$seq),decreasing=TRUE)
tab2<-sort(table(i2$seq),decreasing=TRUE)

i1Summary<-data.frame(
  'bar'=names(tab1),
  'count'=as.vector(tab1),
  'inIndex1'=names(tab1) %in% samples$Index1Sequence,
  'inRevCompIndex1'=names(tab1) %in% revComp(samples$Index1Sequence),
  'inIndex2'=names(tab1) %in% samples$Index2Sequence,
  'inRevCompIndex2'=names(tab1) %in% revComp(samples$Index2Sequence),
stringsAsFactors=FALSE)

i2Summary<-data.frame(
  'bar'=names(tab2),
  'count'=as.vector(tab2),
  'inIndex1'=names(tab2) %in% samples$Index1Sequence,
  'inRevCompIndex1'=names(tab2) %in% revComp(samples$Index1Sequence),
  'inIndex2'=names(tab2) %in% samples$Index2Sequence,
  'inRevCompIndex2'=names(tab2) %in% revComp(samples$Index2Sequence),
stringsAsFactors=FALSE)

write.csv(i1Summary[i1Summary$count>2,],'i1.csv')
write.csv(i2Summary[i2Summary$count>2,],'i2.csv')
