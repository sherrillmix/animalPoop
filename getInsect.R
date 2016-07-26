library(parallel)
library(dnar)
library(levenR)
source('functions.R')

setwd('data/bug')
sras<-c( "SRR611639", "SRR616202", "SRR616203", "SRR616204", "SRR616205", "SRR616206")
for(ii in sras){
	if(!file.exists(sprintf('%s.fastq.gz',ii))){
		cmd<-sprintf('wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/%s/%s/%s.sra',substring(ii,1,6),ii,ii)
		message(cmd)
		system(cmd)
		cmd<-sprintf('fastq-dump -gzip %s.sra',ii)
		message(cmd)
		system(cmd)
	}
}
setwd('../..')


bugs<-read.csv('data/bug/bugs.csv',stringsAsFactors=FALSE)
bugs$sample<-sub(' +$','',bugs$sample)
bugs$sample<-toupper(bugs$sample)
rownames(bugs)<-bugs$sample

runs<-list.files('data/bug','run[1-4].csv',full.names=TRUE)
info<-lapply(runs,read.csv,stringsAsFactors=FALSE)
nSamples<-sapply(info,nrow)
info<-do.call(rbind,info)
info$run<-rep(sub('.*run([0-9]).csv','\\1',runs),nSamples)
#fix various misnamings
info$bak<-info$sample
info$sample<-sub(' +$','',info$sample)
info$sample<-sub('_ ','_',info$sample)
info$sample<-sub(' ','_',info$sample)
info$sample<-toupper(info$sample)
info[info$sample=='DS_2_49','sample']<-'DS_2_049'
info[info$sample=='CH_11','sample']<-'CH_011'
info[info$sample=='CH_11','sample']<-'CH_011'
info[info$sample=='J_191_1','sample']<-'J_191'
info[info$sample=='3_J_204_1_1','sample']<-'J_204_1'
info[info$sample=='11_DS_F_016_1','sample']<-'DS_F_016'
info[!info$sample %in% bugs$sample&grepl('^[0-9]+_',info$sample),'sample']<-sprintf('%s',sub('^[0-9]+_','',info[!info$sample %in% bugs$sample&grepl('^[0-9]+_',info$sample),'sample']))
#no information. just deleting
info<-info[info$sample!='J_331',]

tmp<-info[info$run=='4',]
tmp$run<-'5'
info<-rbind(info,tmp)
info$sample_run<-sprintf('%s__%s',info$sample,info$run)

if(!exists('allSeqs')){
	allSeqs<-do.call(rbind,mclapply(list.files('data/bug','fastq.gz$'),function(x){
		seqs<-read.fastq(sprintf('data/bug/%s',x),convert=TRUE)
		seqs$file<-x
		return(seqs)
	},mc.cores=5))
	allSeqs$bar<-substring(allSeqs$seq,5,12)
	allSeqs$link<-substring(allSeqs$seq,13,14)
	allSeqs$primer<-substring(allSeqs$seq,15,34)
	allSeqs$key<-substring(allSeqs$seq,1,4)
}


#it looks like the files goes in order with 616204 and 616205 duplicated?
fileOrder<-sort(list.files('data/bug','fastq.gz$'))
info$file<-fileOrder[as.numeric(info$run)]

dir.create('work/data/bug',showWarnings=FALSE)
info$name<-info$sample_run
info$species<-sub('[ \n]*\\(.*$','',bugs[info$sample,'name'])
info$order<-bugs[info$sample,'order']
info$diet<-bugs[info$sample,'diet']
info$habitat<-bugs[info$sample,'habitat']
info$growth<-bugs[info$sample,'growth']
info$weight<-sapply(strsplit(bugs[info$sample,'length'],'[~-]'),function(x)mean(as.numeric(x))) #CAREFUL THIS IS LENGTH NOT WEIGHT
info<-info[info$growth=='adult',]

write.csv(info[,c('name','species','diet','habitat','weight')],'work/data/bug/info.csv')

#not sure why TCAG on last run but looks like not a match anyway
seqs<-allSeqs[allSeqs$key=='GACT',]
seqs<-seqs[seqs$link %in% c(info$linkerF,info$linkerR),]
seqs$forward<-seqs$link %in% info$linkerF
seqs<-seqs[seqs$bar %in% c(info$barcodeF,info$barcodeR),]
cbind(table(info$barcodeF,info$run),table(seqs$bar,seqs$file))
seqs<-seqs[paste(seqs$file,seqs$bar) %in% paste(info$file,info$barcodeF),]
seqs$primer<-substring(seqs$seq,15,ifelse(seqs$forward,34,33))
seqs<-seqs[seqs$primer %in% c(info$primerFSeq,info$primerRSeq),]
#lots of reads end in N
seqs$trim<-sub('N.*$','',substring(seqs$seq,ifelse(seqs$forward,34,33)+1))
seqs<-seqs[nchar(seqs$trim)>200&nchar(seqs$trim<600),]

rownames(info)<-paste(info$barcodeF,info$file)
seqs$sample<-info[paste(seqs$bar,seqs$file),'sample_run']

tmp<-runSwarm(unlist(seqs$trim),swarmBin='~/installs/swarm/swarm',swarmArgs='-f -t 32')
otus<-tmp[['otus']]
otuSeqs<-tmp[['seqs']]
write.fa(1:length(otuSeqs),otuSeqs,'work/data/bug/swarmSeqs.fa.gz')
otuTab<-table(seqs$sample,otus)
write.csv(otuTab,'work/data/bug/otuTab.csv')


