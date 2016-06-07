targets<-readLines('data/muegge/desired.txt')
fastqs<-list.files('data/muegge','fastq.gz',full.names=TRUE)
firstLines<-sapply(fastqs,readLines,n=1)
matches<-lapply(targets,grep,firstLines)
if(any(sapply(matches,length)>1))stop(simpleError('Multiple matches found'))
out<-data.frame('name'=targets,'fastq'=fastqs[unlist(matches)])
write.csv(out,'targetFiles.csv',row.names=FALSE)

