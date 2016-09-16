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
logAxis<-function(axisNum=2,exponent=TRUE,addExtra=TRUE,minorTcl=-.2,axisMin=-Inf,...){
  if(axisNum %in% c(2,4)){
    minX<-max(10^par('usr')[3],axisMin)
    maxX<-10^par('usr')[4]
  }else{
    minX<-max(10^par('usr')[1],axisMin)
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
qiimeTaxaToTable<-function(taxa,targetRanks=c("k","p","c","o","f","g","s")){
  unassigned<-taxa=='Unassigned'
  splits<-strsplit(taxa,'; ')
  ranks<-lapply(splits,function(x){
    labels<-sub('^.__','',x)
    names(labels)<-sub('__.*','',x)
    return(labels)
  })
  if(any(nchar(unlist(lapply(ranks[!unassigned],names)))>1))stop(simpleError('Unexpected taxonomic rank found'))
  out<-do.call(rbind,lapply(ranks,function(x)x[targetRanks]))
  colnames(out)<-targetRanks
  out[out=='']<-NA
  return(out)
}
taxaTableToName<-function(taxa,addGenusToSpecies=TRUE,makeUnique=TRUE){
  if(addGenusToSpecies)taxa[,'s']<-ifelse(is.na(taxa[,'s']),taxa[,'s'],paste(taxa[,'g'],taxa[,'s']))
  out<-apply(taxa,1,function(x)paste(colnames(taxa),x,sep='__')[max(c(1,which(!is.na(x))))])
  if(makeUnique)out<-paste(out,ave(out,out,FUN=function(x)1:length(x)))
  return(out)
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
  qiimeTaxaDf<-mclapply(paste(dataDir,studies,sep='/','qiimeTaxa.csv'),read.csv,row.names=1,stringsAsFactors=FALSE,mc.cores=10)
  qiimeTaxa<-mclapply(qiimeTaxaDf,function(x){
    out<-qiimeTaxaToTable(x$taxa)
    out<-cbind(out,'name'=taxaTableToName(out))
    rownames(out)<-rownames(x)
    return(out)
  },mc.cores=10)
  names(qiimeTaxa)<-names(qiimeTaxaDf)<-studies

  if(any(mapply(function(x,y)any(sort(x)!=sort(y)),lapply(qiimeOtus,colnames),lapply(qiimeTaxa,rownames))))stop(simpleError('Mismatch between taxas and OTUs'))
  message('Done reading qiime OTUs')
}
#select otu type here
if(!exists('combinedRare')){
  combinedRare<-lapply(list('qiime'=qiimeOtus,'swarm'=swarmOtus),function(otus){
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
    message('Done Calculating rarefaction')
    return(list('rares'=rares,'minRareAll'=minRareAll))
  })
}

if(!exists('otuProps')){
  otuProps<-lapply(names(qiimeOtus),function(study){
    x<-qiimeOtus[[study]]
    rownames(x)<-paste(info[rownames(x),'species'],ave(1:nrow(x),info[rownames(x),'species'],FUN=function(z)1:length(z)))
    colnames(x)<-qiimeTaxa[[study]][colnames(x),'name']
    out<-t(apply(x,1,function(x)x/sum(x)))
    return(out)
  })
  names(otuProps)<-names(qiimeOtus)

  dummy<-mclapply(unique(info$study),function(ii){
    message(ii)
    thisOtus<-otuProps[[ii]][,apply(otuProps[[ii]],2,max)>.01]
    thisOtus<-as.matrix(log(thisOtus+1))
    pdf(sprintf('out/heat/qiime_%s.pdf',ii),height=20,width=20)
      heatmap(thisOtus,col=rev(heat.colors(100)),main=ii,scale='none',margins=c(8,5))
    dev.off()
    message('Done ',ii)
  },mc.cores=10)

}

for(grouper in c('swarm','qiime')){
  message(grouper)
  if(grouper=='qiime') otus<-qiimeOtus
  else if(grouper=='swarm') otus<-swarmOtus

  rares<-combinedRare[[grouper]][['rares']]
  minRareAll<-combinedRare[[grouper]][['minRareAll']]
  rareAll<-minRareAll['1',]

  info$rare<-rares[info$name]
  info$rareAll<-rareAll[info$name]
  info$aveRare<-ave(info$rare,info$species,info$study,FUN=function(x)mean(x,na.rm=TRUE))
  info$aveRareAll<-ave(info$rareAll,info$species,info$study,FUN=function(x)mean(x,na.rm=TRUE))
  pdf(sprintf('out/%s_raw.pdf',grouper),width=8,height=9)
  par(mfrow=c(3,4))
  for(ii in unique(info$study)){
    thisData<-info[info$study==ii,]
    plot(thisData$weight,thisData$rare,log='xy',main=ii,xlab='Weight (kg)',ylab='Rarefied species',xaxt='n')
    logAxis(1,minorTcl=0,addExtra=FALSE)
    text(thisData$weight,thisData$rare,thisData$species,cex=.5,col='#00000044')
  }
  dev.off()
  pdf(sprintf('out/%s_mean.pdf',grouper),width=8,height=9)
  par(mfrow=c(3,4))
  xlim<-range(info$weight)
  for(ii in c(unique(info$study),'all')){
    if(ii=='all'){
      thisData<-unique(info[,c('weight','aveRareAll','species')])
      colnames(thisData)[colnames(thisData)=='aveRareAll']<-'aveRare'
    } else {
      thisData<-unique(info[info$study==ii,c('weight','aveRare','species')])
    }
    plot(thisData$weight,thisData$aveRare,log='xy',xlab='Weight (kg)',ylab='Rarefied species',xlim=xlim,xaxt='n',las=1)
    logAxis(1,minorTcl=0)
    text(thisData$weight,thisData$aveRare,thisData$species,cex=.5,col='#00000044')
    mod<-addLm(log(thisData$aveRare),log(thisData$weight))
    p<-coef(summary(mod))['lWeight','Pr(>|t|)']
    beta<-coef(mod)['lWeight']
    title(main=sprintf('%s\nr2=%.02f p=%.03f B=%.02f',ii,summary(mod)$r.squared,p,beta))
  }
  dev.off()

  pdf(sprintf('out/%s_minObs.pdf',grouper),width=4,height=4)
    for(ii in as.numeric(rownames(minRareAll))){
      thisData<-info[,c('weight','species')]
      thisData$aveRare<-ave(minRareAll[as.character(ii),],info$species)
      thisData<-unique(thisData)
      par(mar=c(3.1,3,1,.5))
      plot(thisData$weight,thisData$aveRare,log='xy',xlab='Weight (kg)',ylab='Rarefied species',xlim=xlim,las=1,xaxt='n',mgp=c(2,.7,0),yaxt='n')
      logAxis(1,mgp=c(3,.7,0))
      logAxis(2,mgp=c(3,.7,0),exponent=FALSE,las=1)
      text(thisData$weight,thisData$aveRare,thisData$species,cex=.3,col='#00000044')
      mod<-addLm(log(thisData$aveRare),log(thisData$weight))
      title(main=sprintf('>=%d reads',ii))
    }
  dev.off()

}

cols<-rainbow.lab(length(unique(info$study)),alpha=.7)
names(cols)<-unique(info$study)
pdf('out/qiime_vs_swarm.pdf',width=4,height=4)
  qMinRareAll<-combinedRare[['qiime']][['minRareAll']]
  sMinRareAll<-combinedRare[['swarm']][['minRareAll']]
  plot(qMinRareAll['1',],sMinRareAll['1',],xlab='uclust rarefied species count',ylab='swarm rarefied species count',bg=cols[info$study],pch=21,cex=.6,col=NA,las=1)
  abline(0,1,col='#FF000033',lty=2)
  legend('bottomright',names(cols),pch=21,pt.bg=cols,cex=.55,inset=.01,col=NA)
  plot(qMinRareAll['2',],sMinRareAll['2',],xlab='uclust rarefied species count',ylab='swarm rarefied species count',main='No singletons',bg=cols[info$study],pch=21,cex=.5,col=NA,las=1)
  abline(0,1,col='#FF000033',lty=2)
dev.off()

nQiime<-lapply(qiimeOtus,apply,2,sum)
nSwarm<-lapply(swarmOtus,apply,2,sum)
pdf('out/qiime_vs_swarmNs.pdf')
  par(mfrow=c(3,4),mar=c(4,4,1,.1),mgp=c(2.5,1,0))
  for(ii in names(nQiime)){
    ylim<-range(unlist(c(nQiime[ii],nSwarm[ii])))
    xlim<-range(sapply(c(list(1),nQiime[ii],nSwarm[ii]),length))
    plot(1,1,type='n',xlim=xlim/1000,ylim=ylim,log='y',las=1,yaxt='n',main=ii,xlab='Cumulative number of OTUs (x1000)',ylab='Read frequency')
    if(ii == names(nQiime)[[1]])legend('topleft',c('uclust','swarm'),lty=1,col=c('#FF000099','#0000FF99'),bty='n')
    logAxis(2,las=1,addExtra=FALSE,axisMin=1)
    x<-1:length(nQiime[[ii]])
    #lines(x/max(x),sort(nQiime[[ii]]),col='#FF000099')
    lines(x/1000,sort(nQiime[[ii]]),col='#FF000099')
    x<-1:length(nSwarm[[ii]]) 
    lines(x/1000,sort(nSwarm[[ii]]),col='#0000FF99')
  }
dev.off()



