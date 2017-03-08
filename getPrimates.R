
info<-read.table('data/primates/map_final_analysis.txt',sep='\t',stringsAsFactors=FALSE,header=TRUE)
seqs<-read.fa('data/primates/merged.fna.gz')
seqs$sample<-sub('_[0-9]+.*$','',seqs$name)
if(any(!seqs$sample %in% rownames(info)))stop(simpleError('Unknown sample found'))
seqs<-seqs[!grepl('[^ATCG]',seqs$seq),]

runOtuForming(substring(seqs$seq,1,min(nchar(seqs$seq))),seqs$sample,'work/data/primates')
