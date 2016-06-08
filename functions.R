runSwarm<-function(seqs,swarmBin='swarm',swarmArgs='-f'){
	if(any(is.na(seqs)))stop(simpleError('NAs in seqs'))
	seqIds<-as.numeric(as.factor(seqs))
	seqCounts<-ave(seqIds,seqIds,FUN=length)
	seqNames<-sprintf('%08d_%d',seqIds,seqCounts)
	readFile<-tempfile()
	outFile<-tempfile()
	seqFile<-tempfile()
	uniqSelector<-!duplicated(seqs)
	write.fa(seqNames[uniqSelector],seqs[uniqSelector],readFile)
	browser()
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
	comboSelector<-selector
	message('Throwing out ',sum(!selector),' sequences containing non ACTG (',sum(selector),' remaining)')
	selector<-nchar(out)>=minLength
	message('Throwing out ',sum(comboSelector&!selector),' sequences shorter than ',minLength,' (',sum(comboSelector&selector),' remaining)')
	comboSelector<-comboSelector&selector
	if(minQual>0&!is.null(quals)){
		selector<-sapply(strsplit(quals,' '),function(x)sum(as.numeric(x)<minQual))<=maxBadQual
		message('Throwing out ',sum(comboSelector&!selector),' poor quality sequences (',sum(selector&comboSelector),' remaining)')
		comboSelector<-comboSelector&selector
	}
	out<-out[comboSelector]
	return(out)
}
