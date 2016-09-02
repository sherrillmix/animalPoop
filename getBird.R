library(dnar)
source('functions.R')

#https://figshare.com/s/fe202600737e11e58a2906ec4bbcf141
#http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4685052/

info<-read.table(textConnection(sub('^#','',readLines('data/bird/cP315.mapping.txt'))),sep='\t',header=TRUE,stringsAsFactors=FALSE)
info$weight<-as.numeric(ifelse(info$Weightg=='nr',NA,info$Weightg))/1000
info$sample<-info$SampleID
info$species<-paste(info$AOU_genus,info$AOU_species)
info$name<-info$SampleID

dir.create('work/data/bird',showWarnings=FALSE)
write.csv(info[!is.na(info$weight),c('sample','species','name','weight')],'work/data/bird/info.csv')

seqs<-read.fa('data/bird/cP315.final.fasta.gz')
seqs$sample<-sub('_.*$','',seqs$name)
if(any(!seqs$sample %in% info$SampleID))stop(simpleError('Unknown sample found'))

allSeqs<-tapply(seqs$seq,seqs$sample,function(xx)filterReads(xx,minLength=70,maxLength=200))

runOtuForming(unlist(allSeqs),rep(info$file,sapply(allSeqs,length)),'work/data/bird')
