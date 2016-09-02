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
weight<-read.csv('data/anteater/weights.csv',stringsAsFactors=FALSE)
rownames(weight)<-weight$species
info$weight<-weight[info$species,'weight']

dir.create('work/data/anteater',showWarnings=FALSE)
info$name<-info$file
write.csv(info[,c('name','species','weight')],'work/data/anteater/info.csv')

#download files
mapply(function(xx,yy)if(!file.exists(yy))download.file(sprintf('http://%s',xx),yy),info$fastq_ftp,info$target)

allSeqs<-mclapply(info$target,function(x){
	tmp<-read.fastq(x,convert=TRUE)
	seqs<-filterReads(tmp$seq,tmp$qual,minQual=10,minLength=median(nchar(tmp$seq)),maxBadQual=3)
	return(seqs)
},mc.cores=16)


runOtuForming(unlist(allSeqs),rep(info$file,sapply(allSeqs,length)),'work/data/anteater')
