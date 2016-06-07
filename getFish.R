#downloaded files from SRA

#sample info from SRA
ncbi<-readLines('data/fish/biosample_result.txt')
sraLines<-grep('SRA:',ncbi)
idLines<-grep('original_id=',ncbi)
hostLines<-grep('/host=',ncbi)
lineDiff<-hostLines-idLines
if(any(lineDiff!=lineDiff[1]))stop(simpleError('Host and id lines not consistent'))

fishFirstLines<-sapply(list.files('data/fish','fastq$',full.names=TRUE),readLines,n=1)
#fishFirstLines<-sub('^@([^.]+)\\.[0-9].*$','\\1',fishFirstLines)
fishFirstLines<-sub('^[^ ]+ ([^_]+)_.*$','\\1',fishFirstLines)


fish<-data.frame(
	'srs'=sub('.*SRA: ','',ncbi[sraLines])
	,'srr'=sub('.*/([^.]+).fastq','\\1',names(fishFirstLines))
	,'species'=sub('.*\\host="([^"]+)"','\\1',ncbi[hostLines])
	,'id'=sub('.*\\original_id="([^"]+)"','\\1',ncbi[idLines])
	,stringsAsFactors=FALSE
)

write.csv(fish,'data/fish/sampleInfo.csv',row.names=FALSE)

