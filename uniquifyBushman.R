##File name: splitBushman.R
##Creation date: Aug 31, 2015
##Last modified: Tue Sep 01, 2015  07:00AM
##Created by: scott
##Summary: Uniqueft big file from Aubrey

#actually this is a mix of data from all kinds of sequencing runs so I'll keep it together and break up after otus

info<-read.csv('data/bushman/combined_mapping.txt')
seqs<-read.fa('data/bushman/seqs.fna.gz',assumeSingleLine=TRUE)
seqs$sample<-sub('_[0-9]+$','',seqs$name)
selector<-!grepl('[^ACTG]',seqs$seq,perl=TRUE)
message('Throwing out ',sum(!selector),' sequences containing non ACTG (',sum(selector),' remaining)')
seqs<-seqs[selector,]
seqFactor<-as.factor(seqs$seq)
seqs$seqId<-as.numeric(seqFactor)
write.csv(seqs[,c('name','sample','seqId')],'data/bushman/seqIds.csv')
seqCounts<-tapply(seqs$seq,seqs$seqId,length)
out<-data.frame('seqId'=as.numeric(names(seqCounts)),'count'=seqCounts)
#out$hash<-sapply(out$seq,digest,'sha1') #a bit slow
out$seq<-levels(seqFactor)[out$seqId]
out$name<-sprintf('%d_%d',out$seqId,out$count)


write.fa(out$name,out$seq,'data/bushman/uniq.fa.gz')


#for(ii in unique(dat$sample)){
	#thisData<-dat[dat$sample==i,]
	#write.fa(thisData$name,thisData$seq,)
#}
