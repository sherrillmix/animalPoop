forward<-read.csv('ForwardBarcodes.csv',header=FALSE,stringsAsFactors=FALSE)
colnames(forward)<-c('name','seq')
rownames(forward)<-forward$name

reverse<-read.csv('ReverseBarcodes.csv',header=FALSE,stringsAsFactors=FALSE)
colnames(reverse)[1:8]<-c('plate','well','name','revCompIlluminaAdaptor','bar','reversePad','linker','primer')
rownames(reverse)<-sprintf('%s__%s',reverse$plate,reverse$well)

samples<-read.csv('Scott_Island_Biogeo_Combined_MFs_V5.csv',stringsAsFactors=FALSE,check.names=FALSE)
colnames(samples)[colnames(samples)=='#SampleID_temp']<-'SampleID_temp'
origCols<-colnames(samples)
samples$fBarName<-sprintf('gc%d',samples$FBarcodeGolayNumber)
samples$fBar<-forward[samples$fBarName,'seq']
samples$rBarName<-sprintf('%s__%s',sub('Plate ','plate',samples$RPrimerPlateNumber),samples$WellPosition)
samples$rBar<-reverse[samples$rBarName,'bar']
if(any(is.na(samples$rBar))|any(is.na(samples$fBar)))stop('Unknown barcode')
if(any(table(paste(samples$fBar,samples$rBar))>1))stop('Duplicate barcode')
if(any(table(samples[,'#SampleID'])>1))stop('Duplicate sample ID')
samples$Index1Sequence<-revComp(samples$rBar)
samples$Index2Sequence<-revComp(samples$fBar)

tmp<-textConnection('text','w')
write.csv(samples[,origCols],tmp,row.names=FALSE)
close(tmp)
text[1]<-gsub('"','',text[1])
writeLines(text,'sherrill-mix_islandGut.csv')



