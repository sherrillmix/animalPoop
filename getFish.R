library(dnar)
source("functions.R")
library(parallel)

#downloaded files from SRA


#sample info from SRA
ncbi<-readLines('data/fish/biosample_result.txt')
sraLines<-grep('SRA:',ncbi)
idLines<-grep('original_id=',ncbi)
hostLines<-grep('/host=',ncbi)
lineDiff<-hostLines-idLines
if(any(lineDiff!=lineDiff[1]))stop(simpleError('Host and id lines not consistent'))

fishFirstLines<-sapply(list.files('data/fish','fastq$',full.names=TRUE),readLines,n=1)
#fishFirstLines<-sub('^@([^.]+)\\.[0-9].*$','\\1',fishFirstLines)
fishFirstLines<-sub('^[^ ]+ ([^_]+)_.*$','\\1',fishFirstLines)


info<-data.frame(
	'srs'=sub('.*SRA: ','',ncbi[sraLines])
	,'srr'=sub('.*/([^.]+).fastq','\\1',names(fishFirstLines))
	,'species'=sub('.*\\host="([^"]+)"','\\1',ncbi[hostLines])
	,'id'=sub('.*\\original_id="([^"]+)"','\\1',ncbi[idLines])
	,stringsAsFactors=FALSE
)

write.csv(info,'data/fish/sampleInfo.csv',row.names=FALSE)
info$file<-sprintf('data/fish/%s.fastq',fish$srr)

allSeqs<-mclapply(info$file,function(x){
	tmp<-read.fastq(x,convert=TRUE)
	seqs<-filterReads(tmp$seq,minLength=300,maxLength=375,minQual=10,maxBadQual=3)
	return(seqs)
},mc.cores=16)

tmp<-runSwarm(unlist(allSeqs),swarmBin='~/installs/swarm/swarm',swarmArgs='-f -t 32')
otus<-tmp[['otus']]
seqs<-tmp[['seqs']]
samples<-rep(info$file,sapply(allSeqs,length))
dir.create('work/data/fish',showWarnings=FALSE)
write.fa(1:length(seqs),seqs,'work/data/fish/swarmSeqs.fa.gz')
otuTab<-table(samples,otus)
write.csv(otuTab,'work/data/fish/otuTab.csv')


