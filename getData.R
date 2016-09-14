
getFiles<-list.files('.','^get.*\\.R')
getFiles<-getFiles[getFiles!='getData.R']
lapply(getFiles,function(ii){
  message(ii)
  #make sure each starts fresh
  source(ii,local=TRUE)
})
