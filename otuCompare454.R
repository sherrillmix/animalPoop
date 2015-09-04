##File name: otuCompare.R
##Creation date: Aug 26, 2015
##Last modified: Mon Aug 31, 2015  07:00AM
##Created by: scott
##Summary: Read in OTUs and compare between species

source("~/scripts/R/dna.R")

targets<-read.csv('data/mgrast/mgrastInfo.csv',stringsAsFactors=FALSE,row.names=1)
targets<-targets[!grepl('CRBM',targets$Metagenome.Name),] #calorie restricted humans
targets<-targets[!grepl('\\.MG$',targets$Metagenome.Name),] #metagenomic I think
targets<-targets[!duplicated(targets$Metagenome.Name),] #whole table seems duplicated once
targets$name<-sub('\\.16S$','',targets$Metagenome.Name)
weights<-read.csv('data/mgrast/mgrastWeights.csv',stringsAsFactors=FALSE) #rough estimates from wikipedia and  http://www.sciencemag.org/content/332/6032/970/suppl/DC1
rownames(weights)<-weights$name
targets<-cbind(targets,weights[targets$name,-1])
targets$otu<-sprintf('data/mgrast/%s.out',targets$MG.RAST.ID)
if(any(!sapply(targets$otu,file.exists)))stop(simpleError('Missing otu file'))
otus<-lapply(lapply(targets$otu,readLines),strsplit,' ')
otuNs<-lapply(otus,function(x)lapply(x,function(y)as.numeric(sapply(strsplit(y,'_'),'[[',2))))
otuSums<-lapply(otuNs,sapply,sum)
names(otuSums)<-targets$name
nReads<-sapply(otuSums,sum)
targets$rareN<-sapply(otuSums,rareEquation,1500)

#pchs<-as.numeric(as.factor(targets$gut))+20
#cols<-rainbow(length(unique(targets$order)),alpha=.8)
#names(cols)<-unique(targets$order)
#bgs<-cols[targets$order]
pdf('out/weight_vs_otu_454.pdf',height=4,width=4)
	par(mar=c(2.8,3.6,.2,.4))
	plot(targets$weight,targets$rareN,log='xy',xlab='',ylab='Number of OTUs (rarefied to 1,500 reads)',las=1,mgp=c(2.7,.5,0),xaxt='n',pch=21,cex=1.2,bg='#00000011',col='#000000CC')#,pch=pchs,bg=bgs)
	axis(2,4:20*100,rep('',17),tcl=-.2)
	title(xlab='Approximate weight (kg)',mgp=c(1.7,1,0))
	axis(1,c(1,10,100,1000,5000),mgp=c(1,.5,0))
	text(targets$weight,targets$rareN*.98,targets$name,cex=.3,col='#00000044',adj=c(.5,1))
dev.off()

summary(step(glm(I(log2(rareN))~I(log2(weight))+diet+gut,data=targets),direction='both'))
summary(glm(I(log2(rareN))~I(log2(weight)),data=targets))
