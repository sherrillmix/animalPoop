library(RCurl)
library(XML)

info<-read.table('data/anteater/ERP003782.txt',stringsAsFactors=FALSE,sep='\t',header=TRUE)
sampleUrl<-'http://www.ebi.ac.uk/ena/data/view/%s&display=xml'
info$species<-sapply(sprintf(sampleUrl,info$sample_accession),function(address)sub('fecal sample from ','',xpathSApply(xmlParse(address),'//DESCRIPTION',xmlValue)))
info$file<-basename(info$fastq_ftp)
write.csv(info[,c('run_accession','species','file')],'data/anteater/info.csv',row.names=FALSE)


#download files
sapply(info$fastq_ftp,function(x)if(!file.exists(sprintf('data/anteater/%s',basename(x))))download.file(sprintf('http://%s',x),sprintf('data/anteater/%s',basename(x))))

