dataDir<-'work/data/'
studies<-list.files(dataDir)
tmp<-lapply(paste(dataDir,studies,sep='/','info.csv'),read.csv,stringsAsFactors=FALSE)
tmp<-mapply(function(x,y){x$study<-y;x},tmp,studies)
info<-do.call(rbind,lapply(tmp,function(x)x[,c('name','species','weight','study')]))
rownames(info)<-info$name

if(!exists('otus')){
	otus<-lapply(paste(dataDir,studies,sep='/','otuTab.csv'),read.csv,row.names=1)
	names(otus)<-studies
	if(any(mapply(function(x,y)any(!y%in%rownames(x)),otus,split(info$name,info$study))))stop(simpleError('Mismatch between info and OTUs'))
}

rares<-mapply(function(otu,targets){
	ns<-apply(otu[targets,],1,sum)
	cut<-floor(quantile(ns,.1)*.9)
	rares<-apply(otu[targets,],1,rareEquation,samples=floor(cut*.5))
	rares[ns<cut]<-NA	
	return(rares)
},otus,split(info$name,info$study))
rares<-structure(unlist(rares),.Names=unlist(lapply(rares,names)))


info$rare<-rares[info$name]
info$aveRare<-ave(info$rare,info$species,info$study,FUN=function(x)mean(x,na.rm=TRUE))
pdf('out/raw.pdf',width=8,height=6)
par(mfrow=c(2,3))
for(ii in unique(info$study)){
	thisData<-info[info$study==ii,]
	plot(thisData$weight,thisData$rare,log='xy',main=ii,xlab='Weight',ylab='Species')
	text(thisData$weight,thisData$rare,thisData$species,cex=.5,col='#00000044')
}
dev.off()


pdf('out/mean.pdf',width=8,height=6)
par(mfrow=c(2,3))
for(ii in unique(info$study)){
	thisData<-unique(info[info$study==ii,c('weight','aveRare','species')])
	plot(thisData$weight,thisData$aveRare,log='xy',main=ii,xlab='Weight',ylab='Species')
	text(thisData$weight,thisData$aveRare,thisData$species,cex=.5,col='#00000044')
}
dev.off()


