runSwarm<-function(seqs,swarmBin='swarm',swarmArgs='-f'){
	seqIds<-as.numeric(as.factor(seqs))
	seqCounts<-ave(seqIds,seqIds,FUN=length)
	seqNames<-sprintf('%08d_%d',seqIds,seqCounts)
	readFile<-tempfile()
	outFile<-tempfile()
	uniqSelector<-!duplicated(seqs)
	write.fa(seqNames[uniqSelector],seqs[uniqSelector],readFile)
	cmd<-sprintf('%s %s %s -o %s',swarmBin,swarmArgs,readFile,outFile)
	system(cmd)
	swarm<-readLines(outFile)
	otuSplit<-strsplit(swarm,' ')
	otus<-rep(1:length(otuSplit),sapply(otuSplit,length))
	names(otus)<-unlist(otuSplit)	
	out<-otus[seqNames]
	return(out)
}
