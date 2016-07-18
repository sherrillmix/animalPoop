library(dnar)
source("functions.R")
library(parallel)

#downloaded files from SRA

desiredSamples<-c('Gazelle2','Okapi1','BlackRhino1','BigHornW3','Kroo1','SpgbkW','Rabbit','HyraxSTL','BaboonW','Kroo3','Zebra2','Gazelle3','Armadillo','BigHornSD','Wildass1','AfElphSD3','Okapi3','ZebraSTL1','Callimicos','BlackRhino2','Orang1','Okapi2')

weights<-read.csv('data/muegge/targetWeights.csv')
if(any(!desiredSamples %in% weights$name))stop(simpleError('Missing sample'))
info<-read.csv('data/muegge/sra_result.csv',stringsAsFactors=FALSE)
info$name<-sapply(strsplit(info$Library.Name,':'),'[[',1)
if(any(!desiredSamples %in% info$name))stop(simpleError('Missing sample'))
info<-info[info$name %in% desiredSamples,]


firstLines<-sapply(list.files('data/muegge','fastq.gz$',full.names=TRUE),readLines,n=1)
#fishFirstLines<-sub('^@([^.]+)\\.[0-9].*$','\\1',fishFirstLines)
firstLines<-sub('^[^ ]+ ([^.]+)\\..*$','\\1',firstLines)
if(any(!desiredSamples %in% firstLines))stop(simpleError('Missing sample'))
info$file<-sapply(info$name,function(x)names(firstLines)[firstLines==x])
info<-merge(info[,c('name','file')],weights)
dir.create('work/data/mueggeIllumina',showWarnings=FALSE)
info$common<-info$name
info$name<-info$file
write.csv(info,'work/data/mueggeIllumina/info.csv')

allSeqs<-mclapply(info$file,function(x){
	tmp<-read.fastq(x,convert=TRUE)
	seqs<-filterReads(tmp$seq,tmp$qual,minLength=100,minQual=20,maxBadQual=3)
	return(seqs)
},mc.cores=16)

tmp<-runSwarm(unlist(allSeqs),swarmBin='~/installs/swarm/swarm',swarmArgs='-f -t 32')
otus<-tmp[['otus']]
seqs<-tmp[['seqs']]
samples<-rep(info$file,sapply(allSeqs,length))
write.fa(1:length(seqs),seqs,'work/data/mueggeIllumina/swarmSeqs.fa.gz')
otuTab<-table(samples,otus)
write.csv(otuTab,'work/data/mueggeIllumina/otuTab.csv')



