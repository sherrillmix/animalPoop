
samples<-read.csv('sherrill-mix_islandGut.csv',stringsAsFactors=FALSE)
if(any(table(samples$X.SampleID)>1))stop("Duplicate sample")
if(any(table(paste(samples$Index1Sequence,samples$Index2Sequence))>1))stop("Duplicate barcodes")
writeLines(apply(samples[,c('X.SampleID','Index1Sequence','Index2Sequence')],1,paste,collapse=','),'bars.csv')

system('splitbarcodes data/Undetermined_S0_L001_R1_001.fastq.gz data/Undetermined_S0_L001_R2_001.fastq.gz -i data/Undetermined_S0_L001_I1_001.fastq.gz data/Undetermined_S0_L001_I2_001.fastq.gz -o split -b bars.csv -d 10000')


