runSwarm<-function(seqs,swarmBin='swarm',swarmArgs='-f'){
	seqIds<-as.numeric(as.factor(seqs))
	seqCounts<-ave(seqIds,seqIds,FUN=length)
	seqNames<-sprintf('%08d_%d',seqIds,seqCounts)
	readFile<-tempfile()
	outFile<-tempfile()
	seqFile<-tempfile()
	uniqSelector<-!duplicated(seqs)
	write.fa(seqNames[uniqSelector],seqs[uniqSelector],readFile)
	cmd<-sprintf('%s %s %s -o %s -w %s',swarmBin,swarmArgs,readFile,outFile,seqFile)
	system(cmd)
	swarm<-readLines(outFile)
	otuSplit<-strsplit(swarm,' ')
	otus<-rep(1:length(otuSplit),sapply(otuSplit,length))
	names(otus)<-unlist(otuSplit)	
	out<-otus[seqNames]
	seedSeqs<-read.fa(seqFile)
	return(list('otus'=out,'seqs'=seedSeqs))
}

filterReads<-function(seqs,quals=NULL,minLength=0,maxLength=1e6,minQual=0,maxBadQual=3){
	if(maxLength==Inf)maxLength<-max(nchar(seqs)) #substring doesn't like Inf
	out<-sub('[actgn]+$','',sub('^[actgn]+','',seqs)) #remove trimmed sequence at start and end
	out<-substring(out,1,maxLength)
	selector<-!grepl('[^ACTG]',out)
	message('Throwing out ',sum(!selector),' sequences containing non ACTG (',sum(selector),' remaining)')
	out<-out[selector]
	selector<-nchar(out)>=minLength
	message('Throwing out ',sum(!selector),' sequences shorter than ',minLength,' (',sum(selector),' remaining)')
	out<-out[selector]
	if(minQual>0&!is.null(quals)){
		selector<-sapply(strsplit(quals,' '),function(x)sum(as.numeric(x)<minQual))<=maxBadQual
		message('Throwing out ',sum(!selector),' poor quality sequences (',sum(selector),' remaining)')
		out<-out[selector]
	}
	return(out)
}
