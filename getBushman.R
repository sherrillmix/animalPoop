library(dnar)

seqs<-read.fa('data/bushman/seqs.fna.gz',assumeSingleLine=TRUE)
seqs$sample<-sub('_[0-9]+$','',seqs$name)
for(ii in unique(seqs$sample)){
	message(ii)
	write.fa(seqs$name,seqs$seq,sprintf('data/bushman/split/%s.fa.gz',ii))
}

