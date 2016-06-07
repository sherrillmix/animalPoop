#downloaded files from SRA

#sample info from SRA
ncbi<-readLines('data/fish/biosample_result.txt')
sraLines<-grep('SRA:',ncbi)
hostLines<-grep('/host=',ncbi)
lineDiff<-hostLines-sraLines
if(any(lineDiff!=lineDiff[1]))stop(simpleError('Host and sra lines not consistent'))
fish<-data.frame('sra'=sub('.*SRA: ','',ncbi[sraLines]),'species'=sub('.*\\host=','',ncbi[hostLines]),stringsAsFactors=FALSE)
write.csv(fish,'data/fish/sampleInfo.csv',row.names=FALSE)

