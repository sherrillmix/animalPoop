library(dnar)

forward<-read.csv('ForwardBarcodes.csv',header=FALSE,stringsAsFactors=FALSE)
colnames(forward)<-c('name','seq')
rownames(forward)<-forward$name

reverse<-read.csv('ReverseBarcodes2.csv',header=FALSE,stringsAsFactors=FALSE)
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
writeLines(text,'sherrill-MixLauder_islandGut_run1.csv')

writeLines(samples[grep('Chimpanzee|Human|Macaque|Mandrill|Colobus|Gorilla|Bonobo|Mangabey|Primate',samples[,'Animal name']),'#SampleID'],'work/primates.txt')

samples2<-read.csv('sherrill-MixLauder_islandGut2-APL4.csv',stringsAsFactors=FALSE,check.names=FALSE)
samples2<-samples2[!is.na(samples2$sample_name)&samples2$sample_name!='',]
origCols<-colnames(samples2)
samples2$fBarName<-sprintf('gc%d',samples2[,'forward_barcode_i5_index2'])
samples2$fBar<-forward[samples2$fBarName,'seq']
if(any(samples2$fBar!=samples2$i2_sequence))stop('Manually assigned forward barcodes do not match name')
samples2$rBarName<-sprintf('%s__%s',sub('Plate ','plate',samples2$reverse_primer_plate_number),samples2$reverse_barcode_well_i7_index1)
samples2$rBar<-reverse[samples2$rBarName,'bar']
if(any(is.na(samples2$rBar))|any(samples$rBar=='')|any(is.na(samples2$fBar))|any(samples$fBar==''))stop('Unknown barcode')
if(any(table(paste(samples2$fBar,samples2$rBar))>1))stop('Duplicate barcode')
if(any(table(samples2[,'sample_name'])>1))stop('Duplicate sample ID')
samples2$i1_sequence<-samples2$rBar
samples2$barcode_sequence<-paste(samples2$i2_sequence,samples2$i1_sequence,sep='')

write.csv(samples2[,origCols],'sherrill-MixLauder_islandGut_run2.csv',row.names=FALSE)




