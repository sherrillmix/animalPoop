##File name: otuCompare.R
##Creation date: Aug 26, 2015
##Last modified: Tue Sep 01, 2015  10:00AM
##Created by: scott
##Summary: Read in OTUs and compare between species

source("~/scripts/R/dna.R")

targets<-read.table('data/bushman/combined_mapping.txt',sep='\t',stringsAsFactors=FALSE)
colnames(targets)[c(1,5)]<-c('sample','species')
reads<-read.csv('data/bushman/seqIds.csv',row.names=1,stringsAsFactors=FALSE)
otus<-strsplit(readLines('data/bushman/uniq.out'),' ')
otus<-lapply(otus,function(x)sapply(strsplit(x,'_'),function(x)as.numeric(x[1])))
otuLookup<-rep(1:length(otus),sapply(otus,length))
names(otuLookup)<-unlist(otus)
reads$otu<-otuLookup[as.character(reads$seqId)]
otuTable<-table(reads$sample,reads$otu)
targets$nReads<-apply(otuTable[targets$sample,],1,sum)
targets<-targets[targets$nReads>10000,]
targets$rareN<-apply(otuTable[targets$sample,],1,rareEquation,5000)
table(targets$species)

targets[grepl('human baby',targets$species),'species']<-'baby'

animalWeights<-c(
	'drosophila'=.25/1e6, #http://www.genetics.org/content/53/2/237.full.pdf
	'mouse'=18/1e3,
	'rat'=.4, #http://www.criver.com/files/pdfs/rms/us-model-pricing/rm_rm_c_wistar_rats.aspx
	'fish'=1,
	#'breastfed human baby'=10, #pulling 10 out of butt might be better to use adult or exclude #using ~ non developed adult human average #http://www.biomedcentral.com/content/pdf/1471-2458-12-439.pdf
	#'formula a human baby'=10,
	#'formula b human baby'=10,
	'baby'=10,
	'human'=60,
	'HUMPBACK'=36000,
	'FIN'=60000, #http://www.nmfs.noaa.gov/pr/species/mammals/cetaceans/finwhale.htm
	'RIGHT'=55000
)
targets$weight<-animalWeights[targets$species]




speciesAverages<-tapply(targets$rareN,targets$species,function(x)exp(mean(log(x))))


#pchs<-as.numeric(as.factor(targets$gut))+20
#cols<-rainbow(length(unique(targets$order)),alpha=.8)
#names(cols)<-unique(targets$order)
#bgs<-cols[targets$order]
pdf('out/weight_vs_otu_bushman.pdf',height=5,width=6)
	par(mar=c(2.8,3.6,3.2,.4))
	plot(targets$weight,targets$rareN,log='xy',xlab='',ylab='Number of OTUs (rarefied to 5,000 reads)',las=1,mgp=c(2.7,.5,0),xaxt='n',pch=21,cex=1.2,bg='#00000011',col='#000000CC')#,pch=pchs,bg=bgs)
	axis(2,4:20*100,rep('',17),tcl=-.2)
	title(xlab='Approximate weight (kg)',mgp=c(1.7,1,0))
	axis(1,c(.000001,.0001,.01,1,100,10000),c('.000001','.0001','.01','1','100','10000'),mgp=c(1,.5,0))
	axis(3,animalWeights,tolower(names(animalWeights)),las=2,cex=.5,mgp=c(3,.5,0),tcl=-.4)
	text(targets$weight,targets$rareN*.98,targets$name,cex=.3,col='#00000044',adj=c(.5,1))
	segments(animalWeights[names(speciesAverages)]*.5,speciesAverages,animalWeights[names(speciesAverages)]*2,speciesAverages,pch=21,col='#FF000088',lwd=3)
dev.off()
