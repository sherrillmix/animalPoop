##File name: getMgrastData.R
##Creation date: Aug 27, 2015
##Last modified: Mon Jul 11, 2016  10:00AM
##Created by: scott
##Summary: Download Whale data from MGRAST. More difficult that it should be

library(RCurl)
library(XML)
library(parallel)
source('functions.R')
library(dnar)

projectIds<-3854
#The awesome mg-rast javascript in the table so we cant scrape it. saving the pages manually
#ftp://ftp.metagenomics.anl.gov/data/manual/mg-rast-manual.pdf
#ftp://ftp.metagenomics.anl.gov/projects/README.ftp
#http://metagenomics.anl.gov/metagenomics.cgi?page=MetagenomeProject&project=116
#projectUrls<-sprintf('http://metagenomics.anl.gov/metagenomics.cgi?page=MetagenomeProject&project=%d',projectIds)

if(exists('info'))rm(info)
for(ii in projectIds){
	message('Project ',ii)
	html<-htmlParse(sprintf('data/mgrast/mgrastPages/%d.html',ii))
	links<-unlist(xpathSApply(html,'//a/@href'))
	samples<-unique(sub('.*=([0-9.]+)$','\\1',links[grepl('metagenome=[0-9.]+$',links)]))
	sampleUrls<-sprintf('http://api.metagenomics.anl.gov//download/mgm%s?file=050.1',samples)
	names(sampleUrls)<-samples
	thisInfo<-read.table(sprintf('data/mgrast/mgrastPages/%d.tsv',ii),sep='\t',header=TRUE,stringsAsFactors=FALSE)
	for(jj in samples){
		message('  Sample ',jj)
		if(thisInfo[thisInfo$MG.RAST.ID==jj,'Sequence.Type']=='WGS'){
			message('    ',jj,' is WGS. Skipping')
			next()
		}
		#outPath<-sprintf('data/mgrast/%s.fa.gz',jj)
		outPath<-sprintf('data/mgrast/%s.fastq.gz',jj)
		if(!file.exists(outPath)){
			message('    Downloading file')
			fasta<-readLines(sampleUrls[jj])
			outFile<-gzfile(outPath,'w')
			writeLines(fasta,outFile)
			close(outFile)
		}else{
			message('    File exists')
		}
	}
	thisInfo<-thisInfo[thisInfo$Sequence.Type!='WGS',]
	if(!exists('info')){
		info<-thisInfo
	} else {
		sharedCols<-intersect(colnames(info),colnames(thisInfo))
		info<-rbind(info[,sharedCols],thisInfo[,sharedCols])
	}
}
write.csv(info,'data/mgrast/mgrastWhaleInfo.csv')

info$file<-sprintf('data/mgrast/%s.fastq.gz',info$MG.RAST.ID)
allSeqs<-mclapply(info$file,function(x){
	tmp<-read.fastq(x,convert=TRUE)
	seqs<-filterReads(tmp$seq,tmp$qual,minLength=250,maxLength=350,minQual=10,maxBadQual=3)
	return(seqs)
},mc.cores=16)

tmp<-runSwarm(unlist(allSeqs),swarmBin='~/installs/swarm/swarm',swarmArgs='-f -t 32')
otus<-tmp[['otus']]
seqs<-tmp[['seqs']]
samples<-rep(info$file,sapply(allSeqs,length))
dir.create('work/data/whale',showWarnings=FALSE)
write.fa(1:length(seqs),seqs,'work/data/whale/swarmSeqs.fa.gz')
otuTab<-table(samples,otus)
write.csv(otuTab,'work/data/whale/otuTab.csv')
