
samples<-read.csv('sherrill-MixLauder_islandGut.csv',stringsAsFactors=FALSE)
if(any(table(samples$X.SampleID)>1))stop("Duplicate sample")
if(any(table(paste(samples$Index1Sequence,samples$Index2Sequence))>1))stop("Duplicate barcodes")
writeLines(apply(samples[,c('X.SampleID','Index1Sequence','Index2Sequence')],1,paste,collapse=','),'work/bars.csv')

system('splitbarcodes data/Undetermined_S0_L001_R1_001.fastq.gz data/Undetermined_S0_L001_R2_001.fastq.gz -i data/Undetermined_S0_L001_I1_001.fastq.gz data/Undetermined_S0_L001_I2_001.fastq.gz -o split -b work/bars.csv -d 10000')

samples2<-read.csv('sherrill-MixLauder_islandGut2.csv',stringsAsFactors=FALSE)
if(any(table(samples2$sample_name)>1))stop("Duplicate sample")
if(any(table(paste(samples2$i1_sequence,samples2$i2_sequence))>1))stop("Duplicate barcodes")
samples2$i1_sequence<-revComp(samples2$i1_sequence)
samples2$i2_sequence<-revComp(samples2$i2_sequence)
writeLines(apply(samples2[,c('sample_name','i1_sequence','i2_sequence')],1,paste,collapse=','),'work/bars2.csv')

system('splitbarcodes data/run2/Undetermined_S0_L001_R1_001.fastq.gz data/run2/Undetermined_S0_L001_R2_001.fastq.gz -i data/run2/Undetermined_S0_L001_I1_001.fastq.gz data/run2/Undetermined_S0_L001_I2_001.fastq.gz -o split2 -b work/bars2.csv -d 10000')

#library(dnar)
#x1<-read.fastq('data/run2/Undetermined_S0_L001_I1_001.fastq.gz')
#x2<-read.fastq('data/run2/Undetermined_S0_L001_I2_001.fastq.gz')
#y1<-read.fastq('data/Undetermined_S0_L001_I1_001.fastq.gz')
#y2<-read.fastq('data/Undetermined_S0_L001_I2_001.fastq.gz')
#zz<-table(paste(x1$seq,x2$seq,sep='|'))
#x1s<-table(x1$seq)
#x2s<-table(x2$seq)
#y1s<-table(y1$seq)
#y2s<-table(y2$seq)


#all1<-c(samples2[,'i1_sequence'],samples[,'Index1Sequence'])
#all1<-c(all1,revComp(all1))
#all2<-c(samples2[,'i2_sequence'],samples[,'Index2Sequence'])
#all2<-c(all2,revComp(all2))
#allPossible<-c(all1,all2)

#table(x1$seq %in% all1,x2$seq %in% all2)

#table(x1$seq %in% (samples2$i1_sequence),x2$seq %in% (samples2$i2_sequence))
