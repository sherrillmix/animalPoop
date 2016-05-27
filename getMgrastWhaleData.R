##File name: getMgrastData.R
##Creation date: Aug 27, 2015
##Last modified: Tue Sep 29, 2015  07:00AM
##Created by: scott
##Summary: Download Whale data from MGRAST. More difficult that it should be

library(RCurl)
library(XML)

projectIds<-3854
#The awesome mg-rast javascript in the table so we cant scrape it. saving the pages manually
#ftp://ftp.metagenomics.anl.gov/data/manual/mg-rast-manual.pdf
#ftp://ftp.metagenomics.anl.gov/projects/README.ftp
#http://metagenomics.anl.gov/metagenomics.cgi?page=MetagenomeProject&project=116
#projectUrls<-sprintf('http://metagenomics.anl.gov/metagenomics.cgi?page=MetagenomeProject&project=%d',projectIds)

if(exists('info'))rm(info)
for(ii in projectIds){
	message('Project ',ii)
	html<-htmlParse(sprintf('data/mgrastPages/%d.html',ii))
	links<-unlist(xpathSApply(html,'//a/@href'))
	samples<-unique(sub('.*=([0-9.]+)$','\\1',links[grepl('metagenome=[0-9.]+$',links)]))
	sampleUrls<-sprintf('http://api.metagenomics.anl.gov//download/mgm%s?file=050.1',samples)
	names(sampleUrls)<-samples
	thisInfo<-read.table(sprintf('data/mgrastPages/%d.tsv',ii),sep='\t',header=TRUE,stringsAsFactors=FALSE)
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