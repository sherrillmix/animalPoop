library(dnar)

seqs<-read.fa('data/bushman/seqs.fna.gz',assumeSingleLine=TRUE)
seqs$sample<-sub('_[0-9].*','',seqs$name)
for(ii in unique(seqs$sample)){
	message(ii)
	write.fa(seqs[seqs$sample==ii,'name'],seqs[seqs$sample==ii,'seq'],sprintf('data/bushman/split/%s.fa.gz',ii))
}
