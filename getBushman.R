library(dnar)

seqs<-read.fa('data/bushman/raw/seqs.fna.gz',assumeSingleLine=TRUE)
seqs$sample<-sub('_[0-9].*','',seqs$name)

allSeqs<-tapply(seqs$seq,seqs$sample,function(xx)filterReads(xx,minLength=250,maxLength=350))

tmp<-runSwarm(unlist(allSeqs),swarmBin='~/installs/swarm/swarm',swarmArgs='-f -t 32')
otus<-tmp[['otus']]
seqs<-tmp[['seqs']]
samples<-rep(names(allSeqs),sapply(allSeqs,length))
dir.create('work/data/bushman',showWarnings=FALSE)
write.fa(1:length(seqs),seqs,'work/data/bushman/swarmSeqs.fa.gz')
otuTab<-table(samples,otus)
write.csv(otuTab,'work/data/bushman/otuTab.csv')
