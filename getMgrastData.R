library(RCurl)
library(XML)
library(parallel)
source('functions.R')
library(dnar)

projectIds<-113:116
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
	for(jj in samples){
		message('  Sample ',jj)
		outPath<-sprintf('data/mgrast/%s.fa.gz',jj)
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
	thisInfo<-read.table(sprintf('data/mgrast/mgrastPages/%d.tsv',ii),sep='\t',header=TRUE,stringsAsFactors=FALSE)
	if(!exists('info')){
		info<-thisInfo
	} else {
		sharedCols<-intersect(colnames(info),colnames(thisInfo))
		info<-rbind(info[,sharedCols],thisInfo[,sharedCols])
	}
}
info<-info[!duplicated(info$MG.RAST.ID),]
write.csv(info,'data/mgrast/mgrastInfo.csv')

info$file<-sprintf('data/mgrast/%s.fa.gz',info$MG.RAST.ID)
isMeta<-grepl('MG$',info$Metagenome.Name)
info<-info[!isMeta,]

dir.create('work/data/muegge',showWarnings=FALSE)
weights<-read.csv('data/mgrast/mgrastWeights.csv',stringsAsFactors=FALSE) #rough estimates from wikipedia and  http://www.sciencemag.org/content/332/6032/970/suppl/DC1
rownames(weights)<-weights$name
info$name<-sub('\\..*$','',info$Metagenome.Name)
info<-cbind(info,weights[info$name,-1])
info$common<-info$name
info$name<-info$file
write.csv(info,'work/data/muegge/info.csv')

allSeqs<-mclapply(info$file,function(x){
	tmp<-read.fa(x)
	seqs<-filterReads(tmp$seq,minLength=150,maxLength=275)
	return(seqs)
},mc.cores=16)

tmp<-runSwarm(unlist(allSeqs),swarmBin='~/installs/swarm/swarm',swarmArgs='-f -t 32')
otus<-tmp[['otus']]
seqs<-tmp[['seqs']]
samples<-rep(info$file,sapply(allSeqs,length))
write.fa(1:length(seqs),seqs,'work/data/muegge/swarmSeqs.fa.gz')
otuTab<-table(samples,otus)
write.csv(otuTab,'work/data/muegge/otuTab.csv')
