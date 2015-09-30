#File name: uniquify.R
##Creation date: Aug 25, 2015
##Last modified: Tue Sep 29, 2015  01:00PM
##Created by: scott
##Summary: Take fasta and return fasta with unique reads and counts

library('digest')
source("~/scripts/R/dna.R")

args<-commandArgs(TRUE)
input<-args[1]
output<-args[2]
minLength<-as.numeric(args[3])
maxLength<-as.numeric(args[4])
if(is.null(input))stop(simpleError('Input undefined'))
if(is.null(output))stop(simpleError('Output undefined'))
if(is.na(minLength))minLength<-0
if(is.na(maxLength))maxLength<-Inf
isFastq<-grepl('\\.fastq',input)

if(isFastq){
	seqs<-read.fastq(input)
}else{
	seqs<-read.fa(input)
}
#selector<-nchar(seqs$seq)>max(nchar(seqs$seq))
#message('Throwing out ',sum(!selector),' short sequences (',sum(selector),' remaining)')
#seqs<-seqs[selector,]
seqs$seq<-sub('[actgn]+$','',sub('^[actgn]+','',seqs$seq)) #remove trimmed sequence at start and end
seqs$seq<-substring(seqs$seq,1,maxLength)
selector<-!grepl('[^ACTG]',seqs$seq)
message('Throwing out ',sum(!selector),' sequences containing non ACTG (',sum(selector),' remaining)')
seqs<-seqs[selector,]
selector<-nchar(seqs$seq)>minLength
#selector<-sapply(strsplit(seqs$qual,' '),function(x)all(as.numeric(x)>20))
message('Throwing out ',sum(!selector),' sequences shorter than ',minLength,' (',sum(selector),' remaining)')
seqs<-seqs[selector,]
seqCounts<-tapply(seqs$seq,seqs$seq,length)
out<-data.frame('seq'=names(seqCounts),'count'=seqCounts,stringsAsFactors=FALSE)
out$hash<-sapply(out$seq,digest,'sha1')
out$name<-sprintf('%s_%d',substring(out$hash,1,20),out$count)

write.fa(out$name,out$seq,output)




