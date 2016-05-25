library("dnar")

args<-commandArgs(TRUE)
input<-args[1]
output<-args[2]
minLength<-as.numeric(args[3])
maxLength<-as.numeric(args[4])
qualityCut<-as.numeric(args[5])

if(is.null(input))stop(simpleError('Input undefined'))
if(is.null(output))stop(simpleError('Output undefined'))
if(is.na(minLength))minLength<-0
if(is.na(maxLength))maxLength<-Inf
if(is.na(qualityCut))maxLength<-20
isFastq<-grepl('\\.fastq',input)

if(isFastq){
	seqs<-read.fastq(input)
}else{
	seqs<-read.fa(input)
}
seqs$seq<-sub('[actgn]+$','',sub('^[actgn]+','',seqs$seq)) #remove trimmed sequence at start and end
seqs$seq<-substring(seqs$seq,1,maxLength)
selector<-!grepl('[^ACTG]',seqs$seq)
message('Throwing out ',sum(!selector),' sequences containing non ACTG (',sum(selector),' remaining)')
seqs<-seqs[selector,]
selector<-nchar(seqs$seq)>minLength
message('Throwing out ',sum(!selector),' sequences shorter than ',minLength,' (',sum(selector),' remaining)')
seqs<-seqs[selector,]
if(qualityCut>0&isFastq){
	selector<-sapply(strsplit(seqs$qual,' '),function(x)all(as.numeric(x)>qualityCut))
	message('Throwing out ',sum(!selector),' poor quality sequences (',sum(selector),' remaining)')
	seqs<-seqs[selector,]
}
seqCounts<-tapply(seqs$seq,seqs$seq,length)
out<-data.frame('seq'=names(seqCounts),'count'=seqCounts,stringsAsFactors=FALSE)
out$name<-sprintf('%06d_%d',1:nrow(out),out$count)

write.fa(out$name,out$seq,output)

