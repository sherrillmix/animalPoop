library(dnar)
source('functions.R')

info<-read.table('data/bushman/combined_mapping.txt',sep='\t',stringsAsFactors=FALSE)
colnames(info)[c(1,5)]<-c('sample','species')
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
info$weight<-animalWeights[info$species]
info$file<-info$sample
dir.create('work/data/bushman',showWarnings=FALSE)
write.csv(info[,c('sample','file','species','weight')],'work/data/bushman/info.csv')

seqs<-read.fa('data/bushman/raw/seqs.fna.gz',assumeSingleLine=TRUE)
seqs$sample<-sub('_[0-9].*','',seqs$name)

allSeqs<-tapply(seqs$seq,seqs$sample,function(xx)filterReads(xx,minLength=250,maxLength=350))

tmp<-runSwarm(unlist(allSeqs),swarmBin='~/installs/swarm/swarm',swarmArgs='-f -t 32')
otus<-tmp[['otus']]
seqs<-tmp[['seqs']]
samples<-rep(names(allSeqs),sapply(allSeqs,length))
write.fa(1:length(seqs),seqs,'work/data/bushman/swarmSeqs.fa.gz')
otuTab<-table(samples,otus)
write.csv(otuTab,'work/data/bushman/otuTab.csv')
