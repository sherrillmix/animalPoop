##File name: otuCompare.R
##Creation date: Aug 26, 2015
##Last modified: Tue Sep 08, 2015  07:00AM
##Created by: scott
##Summary: Read in OTUs and compare between species

source("~/scripts/R/dna.R")

targets<-read.csv('data/anteater/info.csv',stringsAsFactors=FALSE)
weight<-read.csv('data/anteater/weights.csv',stringsAsFactors=FALSE)
rownames(weight)<-weight$species
targets$weight<-weight[targets$species,'weight']

targets$otu<-sprintf('data/anteater/%s',sub('.fastq.gz$','.out',targets$file))
if(any(!sapply(targets$otu,file.exists)))stop(simpleError('Missing otu file'))
otus<-lapply(lapply(targets$otu,readLines),strsplit,' ')
otuNs<-lapply(otus,function(x)lapply(x,function(y)as.numeric(sapply(strsplit(y,'_'),'[[',2))))
otuSums<-lapply(otuNs,sapply,sum)
names(otuSums)<-targets$run_accession
targets$nReads<-sapply(otuSums,sum)
targets<-targets[targets$nReads>2500,]
targets$rareN<-sapply(otuSums[targets$run_accession],rareEquation,2000)

speciesAverages<-tapply(targets$rareN,targets$species,function(x)exp(mean(log(x))))

pdf('out/weight_vs_otu_anteater.pdf',height=4,width=4)
	par(mar=c(2.8,3.6,.2,.4))
	plot(targets$weight,targets$rareN,log='xy',xlab='',ylab='Number of OTUs (rarefied to 2,000 reads)',las=1,mgp=c(2.7,.5,0),xaxt='n',pch=21,cex=1.2,bg='#00000011',col='#000000CC')
	axis(2,4:20*100,rep('',17),tcl=-.2)
	title(xlab='Approximate weight (kg)',mgp=c(1.7,1,0))
	axis(1,c(1,10,100,1000,5000),mgp=c(1,.5,0))
	text(targets$weight,targets$rareN*.98,targets$species,cex=.35,col='#00000044',adj=c(.5,1))
	segments(weight[names(speciesAverages),'weight']*.8,speciesAverages,weight[names(speciesAverages),'weight']*1.25,speciesAverages,pch=21,col='#FF000088',lwd=3)
dev.off()

