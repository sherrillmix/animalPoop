##File name: findTargets.R
##Creation date: Aug 25, 2015
##Last modified: Wed Aug 26, 2015  08:00AM
##Created by: scott
##Summary: Read in desired names and find in fastq files

targets<-readLines('data/desired.txt')
fastqs<-list.files('data','fastq.gz',full.names=TRUE)
firstLines<-sapply(fastqs,readLines,n=1)
matches<-lapply(targets,grep,firstLines)
if(any(sapply(matches,length)>1))stop(simpleError('Multiple matches found'))
out<-data.frame('name'=targets,'fastq'=fastqs[unlist(matches)])
write.csv(out,'targetFiles.csv',row.names=FALSE)

