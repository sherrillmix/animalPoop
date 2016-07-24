library(parallel)
library(dnar)

setwd('data/bug')
sras<-c( "SRR611639", "SRR616202", "SRR616203", "SRR616204", "SRR616205", "SRR616206")
for(ii in sras){
	if(!file.exists(sprintf('%s.fastq.gz',ii))){
		cmd<-sprintf('wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/%s/%s/%s.sra',substring(ii,1,6),ii,ii)
		message(cmd)
		system(cmd)
		cmd<-sprintf('fastq-dump -gzip %s.sra',ii)
		message(cmd)
		system(cmd)
	}
}
setwd('../..')


bugs<-read.csv('data/bug/bugs.csv')

runs<-list.files('data/bug','run[1-4].csv',full.names=TRUE)
info<-lapply(runs,read.csv,stringsAsFactors=FALSE)
nSamples<-sapply(info,nrow)
info<-do.call(rbind,info)
info$run<-rep(sub('.*run([0-9]).csv','\\1',runs),nSamples)

allSeqs<-do.call(rbind,mclapply(list.files('data/bug','fastq.gz$'),function(x){
	seqs<-read.fastq(sprintf('data/bug/%s',x),convert=TRUE)
	seqs$file<-x
	return(seqs)
},mc.cores=5))
allSeqs$bar<-substring(allSeqs$seq,5,12)
allSeqs$link<-substring(allSeqs$seq,13,14)
allSeqs$primer<-substring(allSeqs$seq,15,34)
allSeqs$key<-substring(allSeqs$seq,1,4)
#not sure why TCAG on last run but looks like not a match anyway
seqs<-allSeqs[allSeqs$key=='GACT',]
seqs<-seqs[seqs$link %in% c(info$linkerF,info$linkerR),]
seqs<-seqs[seqs$bar %in% c(info$barcodeF,info$barcodeR),]


