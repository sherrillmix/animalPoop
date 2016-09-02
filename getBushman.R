library(dnar)
source('functions.R')

info<-read.table('data/bushman/raw/combined_mapping.txt',sep='\t',stringsAsFactors=FALSE)
colnames(info)[c(1,5)]<-c('sample','species')
info<-info[,c('sample','species')]
animalWeights<-c(
  'drosophila'=.25/1e6, #http://www.genetics.org/content/53/2/237.full.pdf
  'mouse'=18/1e3,
  'rat'=.4, #http://www.criver.com/files/pdfs/rms/us-model-pricing/rm_rm_c_wistar_rats.aspx
  'wild rat'=.4,
  'macaque'=6.5, #https://en.wikipedia.org/wiki/Rhesus_macaque
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
dir.create('work/data/bushman',showWarnings=FALSE)
info<-info[!is.na(info$species),]
info<-info[!grepl('baby',info$species),]
info<-info[!grepl('drosophila',info$species)|grepl('w1118',info$sample),]
info<-info[,c('sample','species')]
additional<-read.csv('data/bushman/raw/macaqueRat.csv',stringsAsFactors=FALSE)
additional$species[grepl('wild.[rR]at',additional$sample)]<-'wild rat'
additional$species[additional$species=='Macaca mulatta']<-'macaque'
additional<-additional[additional$sample!='Day14.JK54.feces.JK54',] #not sequenced apparently
info<-rbind(info,additional)

info$weight<-animalWeights[info$species]
info$name<-info$sample
write.csv(info[,c('sample','name','species','weight')],'work/data/bushman/info.csv')

seqs<-readFaDir('data/bushman/raw/','fna.gz$',assumeSingleLine=TRUE)
seqs$sample<-sub('_[0-9].*','',seqs$name)

allSeqs<-tapply(seqs$seq,seqs$sample,function(xx)filterReads(xx,minLength=250,maxLength=350))

runOtuForming(unlist(allSeqs),rep(info$file,sapply(allSeqs,length)),'work/data/bushman')
