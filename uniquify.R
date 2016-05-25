library('digest')
source("~/scripts/R/dna.R")

args<-commandArgs(TRUE)
input<-args[1]
output<-args[2]

qualityCut<-ifelse(length(args)>2,as.numeric(args[3]),20)

if(is.null(input))stop(simpleError('Input undefined'))
if(is.null(output))stop(simpleError('Output undefined'))

seqs<-read.fastq(input)
selector<-nchar(seqs$seq)==max(nchar(seqs$seq))
message('Throwing out ',sum(!selector),' short sequences (',sum(selector),' remaining)')
seqs<-seqs[selector,]
selector<-!grepl('[^ACTG]',seqs$seq)
message('Throwing out ',sum(!selector),' sequences containing non ACTG (',sum(selector),' remaining)')
seqs<-seqs[selector,]
selector<-sapply(strsplit(seqs$qual,' '),function(x)all(as.numeric(x)>qualityCut))
message('Throwing out ',sum(!selector),' poor quality sequences (',sum(selector),' remaining)')
seqs<-seqs[selector,]
seqCounts<-tapply(seqs$seq,seqs$seq,length)
out<-data.frame('seq'=names(seqCounts),'count'=seqCounts)
out$hash<-sapply(out$seq,digest,'sha1')
out$name<-sprintf('%s_%d',substring(out$hash,1,20),out$count)

write.fa(out$name,out$seq,output)




