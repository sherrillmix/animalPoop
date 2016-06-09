library(RCurl)
library(XML)
library(dnar)
library(parallel)
source('functions.R')

info<-read.table('data/anteater/ERP003782.txt',stringsAsFactors=FALSE,sep='\t',header=TRUE)
sampleUrl<-'http://www.ebi.ac.uk/ena/data/view/%s&display=xml'
info$species<-sapply(sprintf(sampleUrl,info$sample_accession),function(address)sub('fecal sample from ','',xpathSApply(xmlParse(address),'//DESCRIPTION',xmlValue)))
info$file<-basename(info$fastq_ftp)
write.csv(info[,c('run_accession','species','file')],'data/anteater/info.csv',row.names=FALSE)
info$target<-sprintf('data/anteater/%s',info$file)


#download files
mapply(function(xx,yy)if(!file.exists(yy))download.file(sprintf('http://%s',xx),yy),info$fastq_ftp,info$target)

allSeqs<-mclapply(info$target,function(x){
	tmp<-read.fastq(x,convert=TRUE)
	seqs<-filterReads(tmp$seq,tmp$qual,minQual=10,minLength=median(nchar(tmp$seq)),maxBadQual=3)
	return(seqs)
},mc.cores=16)

tmp<-runSwarm(unlist(allSeqs),swarmBin='~/installs/swarm/swarm',swarmArgs='-f -t 32')
otus<-tmp[['otus']]
seqs<-tmp[['seqs']]
samples<-rep(info$file,sapply(allSeqs,length))
dir.create('work/data/anteater',showWarnings=FALSE)
write.fa(1:length(seqs),seqs,'work/data/anteater/swarmSeqs.fa.gz')
otuTab<-table(samples,otus)
write.csv(otuTab,'work/data/anteater/otuTab.csv')
