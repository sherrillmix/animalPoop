library(dnar)
library(parallel)

dataDir<-'work/data/'
studies<-list.files(dataDir)
tmp<-lapply(paste(dataDir,studies,sep='/','info.csv'),read.csv,stringsAsFactors=FALSE)
tmp<-mapply(function(x,y){x$study<-y;x},tmp,studies)
info<-do.call(rbind,lapply(tmp,function(x)x[,c('name','species','weight','study')]))
rownames(info)<-info$name

if(!exists('otus')){
	otus<-mclapply(paste(dataDir,studies,sep='/','otuTab.csv'),read.csv,row.names=1,mc.cores=8)
	names(otus)<-studies
	if(any(mapply(function(x,y)any(!y%in%rownames(x)),otus,split(info$name,info$study))))stop(simpleError('Mismatch between info and OTUs'))
}

rares<-mapply(function(otu,targets){
	ns<-apply(otu[targets,],1,sum)
	cut<-max(200,floor(quantile(ns,.1)*.9))
	rares<-apply(otu[targets,],1,rareEquation,samples=floor(cut*.5))
	rares[ns<cut]<-NA	
	return(rares)
},otus,split(info$name,info$study))
rares<-structure(unlist(rares),.Names=unlist(lapply(rares,names)))


info$rare<-rares[info$name]
info$aveRare<-ave(info$rare,info$species,info$study,FUN=function(x)mean(x,na.rm=TRUE))
pdf('out/raw.pdf',width=8,height=9)
par(mfrow=c(3,3))
for(ii in unique(info$study)){
	thisData<-info[info$study==ii,]
	plot(thisData$weight,thisData$rare,log='xy',main=ii,xlab='Weight',ylab='Species')
	text(thisData$weight,thisData$rare,thisData$species,cex=.5,col='#00000044')
}
dev.off()

pdf('out/mean.pdf',width=8,height=9)
par(mfrow=c(3,3))
xlim<-range(info$weight)
for(ii in unique(info$study)){
	thisData<-unique(info[info$study==ii,c('weight','aveRare','species')])
	#bugs are lengths so try cubing for weight
	if(grepl('bug',ii))thisData$weight<-thisData$weight^3/1e3
	plot(thisData$weight,thisData$aveRare,log='xy',xlab='Weight',ylab='Species',xlim=xlim)
	text(thisData$weight,thisData$aveRare,thisData$species,cex=.5,col='#00000044')
	lRare<-log(thisData$aveRare)
	lWeight<-log(thisData$weight)
	mod<-lm(lRare~lWeight)
	fakeWeights<-exp(-15:15)
	coef<-mod$coefficients
	pred<-exp(predict(mod,data.frame('lWeight'=log(fakeWeights))))
	ci<-predict(mod,data.frame('lWeight'=log(fakeWeights)),interval='conf')[,c('lwr','upr')]
	#pred2<-exp(log(fakeWeights)*coef['lWeight']+coef['(Intercept)'])
	polygon(c(fakeWeights,rev(fakeWeights)),exp(c(ci[,'lwr'],rev(ci[,'upr']))),border=NA,col='#FF000033')
	lines(fakeWeights,pred,col='#FF000088',lty=2)
	p<-coef(summary(mod))['lWeight','Pr(>|t|)']
	beta<-coef(mod)['lWeight']
	title(main=sprintf('%s\nr2=%.02f p=%.03f B=%.02f',ii,summary(mod)$r.squared,p,beta))
}
dev.off()


otuProps<-lapply(otus,function(x)t(apply(x,1,function(x)x/sum(x))))
mclapply(unique(info$study),function(ii){
	message(ii)
	pdf(sprintf('out/heat/%s.pdf',ii),height=20,width=20)
	thisOtus<-otuProps[[ii]][,apply(otuProps[[ii]],2,max)>.01]
	thisOtus<-as.matrix(log(thisOtus+1))
	rownames(thisOtus)<-info[rownames(thisOtus),'species']
	heatmap(thisOtus,col=rev(heat.colors(100)),main=ii,scale='column')
	dev.off()
	message('Done ',ii)
},mc.cores=10)
