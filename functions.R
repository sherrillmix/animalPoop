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

runQiime<-function(seqs){
  if(any(is.na(seqs)))stop(simpleError('NAs in seqs'))
  seqIds<-1:length(seqs)
  readDir<-tempfile()
  dir.create(readDir)
  readFile<-file.path(readDir,'XXX.fa')
  outDir<-tempfile()
  seqNames<-sprintf('XXX_%d',seqIds)
  write.fa(seqNames,seqs,readFile)
  #miniconda doesn't like sh so need to use bash
  cmd<-sprintf('echo "source activate qiime1; pick_de_novo_otus.py --input %s --output %s --parallel --jobs_to_start 16 --force"|bash',readFile,outDir)
  message(cmd)
  exit<-system(cmd)
  if(exit!=0)stop(simpleError('Problem running qiime'))
  #get otu assignments
  assigns<-strsplit(readLines(file.path(outDir,'uclust_picked_otus/XXX_otus.txt')),'\t')
  names(assigns)<-sapply(assigns,'[[',1)
  assigns<-lapply(assigns,'[',-1)
  otus<-rep(names(assigns),sapply(assigns,length))
  names(otus)<-unlist(assigns)
  out<-otus[seqNames]
  #get taxa assignments
  taxa<-strsplit(readLines(file.path(outDir,'uclust_assigned_taxonomy/XXX_rep_set_tax_assignments.txt')),'\t')
  names(taxa)<-sapply(taxa,'[[',1)
  taxa<-sapply(taxa,'[[',2)
  seqs<-read.fa(file.path(outDir,'pynast_aligned_seqs/XXX_rep_set_aligned_pfiltered.fasta'))
  seqs<-structure(seqs$seq,.Names=seqs$name)
  #get sequences
  return(list('otus'=out,'seqs'=seqs,'taxa'=taxa))
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


runOtuForming<-function(seqs,samples,outDir){
  #swarm
  tmp<-runSwarm(seqs,swarmBin='~/installs/swarm/swarm',swarmArgs='-f -t 32')
  otus<-tmp[['otus']]
  otuSeqs<-tmp[['seqs']]
  write.fa(1:length(otuSeqs),otuSeqs,file.path(outDir,'swarmSeqs.fa.gz'))
  otuTab<-table(samples,otus)
  write.csv(otuTab,file.path(outDir,'otuTab.csv'))
  #qiime
  qiime<-runQiime(seqs)
  otus<-qiime[['otus']]
  otuTab<-table(samples,otus)
  write.csv(otuTab,file.path(outDir,'qiimeOtuTab.csv'))
  otuSeqs<-qiime[['seqs']]
  write.fa(names(otuSeqs),otuSeqs,file.path(outDir,'qiimeSeqs.fa.gz'))
  taxa<-qiime[['taxa']]
  write.csv(data.frame('name'=names(taxa),'taxa'=taxa,stringsAsFactors=FALSE),file.path(outDir,'qiimeTaxa.csv'))
}
