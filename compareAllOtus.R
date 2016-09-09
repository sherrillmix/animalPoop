library(dnar)
library(parallel)

dataDir<-'work/data/'
studies<-list.files(dataDir)
tmp<-lapply(paste(dataDir,studies,sep='/','info.csv'),read.csv,stringsAsFactors=FALSE)
tmp<-mapply(function(x,y){x$study<-y;x},tmp,studies)
info<-do.call(rbind,lapply(tmp,function(x)x[,c('name','species','weight','study')]))
rownames(info)<-info$name

addLm<-function(lRare,lWeight){
  mod<-lm(lRare~lWeight)
  fakeWeights<-10^(-15:15)
  coef<-mod$coefficients
  pred<-exp(predict(mod,data.frame('lWeight'=log(fakeWeights))))
  ci<-predict(mod,data.frame('lWeight'=log(fakeWeights)),interval='conf')[,c('lwr','upr')]
  polygon(c(fakeWeights,rev(fakeWeights)),exp(c(ci[,'lwr'],rev(ci[,'upr']))),border=NA,col='#FF000033')
  lines(fakeWeights,pred,col='#FF000088',lty=2)
  return(mod)
}
calcRares<-function(otu,targets,cut=max(200,floor(quantile(ns,.1)*.9)),...){
  ns<-apply(otu[targets,],1,sum)
  rares<-apply(otu[targets,],1,rareEquation,samples=floor(cut*.5),...)
  rares[ns<cut]<-NA  
  return(rares)
}
logAxis<-function(axisNum=2,exponent=TRUE,addExtra=TRUE,minorTcl=-.2,...){
  if(axisNum %in% c(2,4)){
    minX<-10^par('usr')[3]
    maxX<-10^par('usr')[4]
  }else{
    minX<-10^par('usr')[1]
    maxX<-10^par('usr')[2]
  }
  allTicks<-unlist(lapply(floor(log10(minX)):ceiling(log10(maxX)),function(x)1:9*10^x))
  allTicks<-allTicks[allTicks<=maxX & allTicks>=minX]
  axis(axisNum,allTicks,rep('',length(allTicks)),tcl=minorTcl)
  prettyY<-seq(ceiling(log10(minX)),floor(log10(maxX)),1)
  axis(axisNum,10^prettyY,rep('',length(prettyY)),tcl=minorTcl*2)
  if(length(prettyY)>7)prettyY<-pretty(prettyY)
  if(addExtra){
    if(length(prettyY)<5)prettyY<-unique(c(prettyY,prettyY+log10(5),prettyY-log10(2)))
    if(length(prettyY)<5)prettyY<-unique(c(prettyY,prettyY+log10(2),prettyY-log10(5)))
  }
  if(exponent) labs<-ifelse(prettyY==0,1,sapply(prettyY,function(x)as.expression(bquote(10^.(x)))))
  else labs<-10^prettyY
  axis(axisNum,10^prettyY,labs,...)
}


if(!exists('swarmOtus')){
  message('Reading swarm OTUs')
  swarmOtus<-mclapply(paste(dataDir,studies,sep='/','otuTab.csv'),read.csv,row.names=1,mc.cores=8)
  names(swarmOtus)<-studies
  if(any(mapply(function(x,y)any(!y%in%rownames(x)),swarmOtus,split(info$name,info$study))))stop(simpleError('Mismatch between info and OTUs'))
  message('Done reading OTUs')
}
if(!exists('qiimeOtus')){
  message('Reading qiime OTUs')
  qiimeOtus<-mclapply(paste(dataDir,studies,sep='/','qiimeOtuTab.csv'),read.csv,row.names=1,mc.cores=8)
  names(qiimeOtus)<-studies
  if(any(mapply(function(x,y)any(!y%in%rownames(x)),qiimeOtus,split(info$name,info$study))))stop(simpleError('Mismatch between info and OTUs'))
  message('Done reading qiime OTUs')
}
#select otu type here
for(grouper in c('swarm','qiime')){
  message(grouper)
  if(grouper=='qiime') otus<-qiimeOtus
  else if(grouper=='swarm') otus<-swarmOtus

  message('Calculating rarefaction')
  rares<-mapply(calcRares,otus,split(info$name,info$study))
  rares<-structure(unlist(rares),.Names=unlist(lapply(rares,names)))

  nCuts<-1:10
  minRareAll<-do.call(rbind,mclapply(nCuts,function(minObs){
    message(minObs)
    rareAll<-mapply(calcRares,otus,split(info$name,info$study),MoreArgs=list(cut=400,minObs=minObs))
    rareAll<-structure(unlist(rareAll),.Names=unlist(lapply(rareAll,names)))
    return(rareAll)
  },mc.cores=5))
  rownames(minRareAll)<-nCuts
  minRareAll<-minRareAll[,info$name]
  rareAll<-minRareAll['1',]
  message('Done Calculating rarefaction')

  info$rare<-rares[info$name]
  info$rareAll<-rareAll[info$name]
  info$aveRare<-ave(info$rare,info$species,info$study,FUN=function(x)mean(x,na.rm=TRUE))
  info$aveRareAll<-ave(info$rareAll,info$species,info$study,FUN=function(x)mean(x,na.rm=TRUE))
  pdf(sprintf('out/%s_raw.pdf',grouper),width=8,height=9)
  par(mfrow=c(3,3))
  for(ii in unique(info$study)){
    thisData<-info[info$study==ii,]
    plot(thisData$weight,thisData$rare,log='xy',main=ii,xlab='Weight',ylab='Species',xaxt='n')
    logAxis(1,minorTcl=0)
    text(thisData$weight,thisData$rare,thisData$species,cex=.5,col='#00000044')
  }
  dev.off()
  pdf(sprintf('out/%s_mean.pdf',grouper),width=8,height=9)
  par(mfrow=c(3,3))
  xlim<-range(info$weight)
  for(ii in c(unique(info$study),'all')){
    if(ii=='all'){
      thisData<-unique(info[,c('weight','aveRareAll','species')])
      colnames(thisData)[colnames(thisData)=='aveRareAll']<-'aveRare'
    } else {
      thisData<-unique(info[info$study==ii,c('weight','aveRare','species')])
    }
    plot(thisData$weight,thisData$aveRare,log='xy',xlab='Weight',ylab='Species',xlim=xlim,xaxt='n',las=1)
    logAxis(1,minorTcl=0)
    text(thisData$weight,thisData$aveRare,thisData$species,cex=.5,col='#00000044')
    mod<-addLm(log(thisData$aveRare),log(thisData$weight))
    p<-coef(summary(mod))['lWeight','Pr(>|t|)']
    beta<-coef(mod)['lWeight']
    title(main=sprintf('%s\nr2=%.02f p=%.03f B=%.02f',ii,summary(mod)$r.squared,p,beta))
  }
  dev.off()

  pdf(sprintf('out/%s_minObs.pdf',grouper),width=4,height=4)
    for(ii in nCuts){
      thisData<-info[,c('weight','species')]
      thisData$aveRare<-ave(minRareAll[as.character(ii),],info$species)
      thisData<-unique(thisData)
      par(mar=c(3.1,3,1,.5))
      plot(thisData$weight,thisData$aveRare,log='xy',xlab='Weight (kg)',ylab='Rarefied species',xlim=xlim,las=1,xaxt='n',mgp=c(2,.7,0))
      logAxis(1,mgp=c(3,.7,0))
      #logAxis(2,mgp=c(3,.7,0),exponent=FALSE,las=1)
      text(thisData$weight,thisData$aveRare,thisData$species,cex=.3,col='#00000044')
      mod<-addLm(log(thisData$aveRare),log(thisData$weight))
      title(main=sprintf('>=%d reads',ii))
    }
  dev.off()

  otuProps<-lapply(otus,function(x)t(apply(x,1,function(x)x/sum(x))))
  mclapply(unique(info$study),function(ii){
    message(ii)
    pdf(sprintf('out/heat/%s_%s.pdf',grouper,ii),height=20,width=20)
    thisOtus<-otuProps[[ii]][,apply(otuProps[[ii]],2,max)>.01]
    thisOtus<-as.matrix(log(thisOtus+1))
    rownames(thisOtus)<-info[rownames(thisOtus),'species']
    heatmap(thisOtus,col=rev(heat.colors(100)),main=ii,scale='column')
    dev.off()
    message('Done ',ii)
  },mc.cores=10)

}
