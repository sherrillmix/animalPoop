
getFiles<-list.files('.','^get.*\\.R')
getFiles<-getFiles[getFiles!='getData.R']
for(ii in getFiles){
  message(ii)
  source(ii)
  rm(list=ls())
}
