nmfs<-read.csv('shark/SPIRALVALVESAMPLEDATA09-06-16-1.csv')
bush<-read.csv('shark/sharkSample.csv')
colnames(bush)[1:3]<-colnames(nmfs)[1:3]
#inferred from missing specimen ID with same sample and haul
bush[bush$SAMPLE_IDENTIFIER=='ARS523'&bush$HAUL==1&bush$SPECIMEN_NUMBER=='?','SPECIMEN_NUMBER']<-105
merged<-merge(bush,nmfs,all=TRUE)
#print(merged[is.na(merged$DATE),])
if(any(is.na(merged$DATE)))stop(simpleError('Missing sample info'))
print(merged[is.na(merged$species),])
merged<-merged[!is.na(merged$species),]

