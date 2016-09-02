library(dnar)
source("functions.R")
library(parallel)


tmp<-list(
	#http://fishbull.noaa.gov/1003/15nelson.pdf
	#10^(-4.53+3.03*log10(200))
	data.frame("Lagodon rhomboides","Pinfish",276)
	#http://onlinelibrary.wiley.com/doi/10.1002/etc.3239/full
	,data.frame("Fundulus heteroclitus","Mummichog",5.5/1000)
	#https://en.wikipedia.org/wiki/Mahi-mahi
	,data.frame("Coryphaena hippurus",'Mahi-mahi',13)
	#http://sedarweb.org/docs/wsupp/S18-RD46%20NC%20RD%20age,%20growth,%20mortality,%20reproductive%20biol.pdf
	#https://en.wikipedia.org/wiki/Red_drum
	,data.frame("Sciaeops ocellatus",'Red drum',10)
	#https://www.researchgate.net/publication/230273336_Length-weight_relationships_for_24_fish_species_in_a_coastal_lagoon_of_the_Mexican_South_Pacific
	,data.frame("Bairdiella chrysoura",'Silver perch',100/1000)
	#https://en.wikipedia.org/wiki/Great_barracuda
	,data.frame("Sphyraena barracuda",'Great barracuda',9)
	#http://www.mass.gov/eea/agencies/dfg/dmf/recreational-fishing/black-sea-bass.html
	,data.frame("Centropristis striata",'Black sea bass',1)
	#https://books.google.com/books?id=SOlDCgAAQBAJ&pg=PA50&lpg=PA50&dq=Caranx+hippos+weight+average&source=bl&ots=2YDP8bTWzJ&sig=lVajLNGco1C14qFN2k99mXE25Bw&hl=en&sa=X&ved=0ahUKEwizzMLUyPPNAhVFlB4KHeaCDu4Q6AEISDAG#v=onepage&q=Caranx%20hippos%20weight%20average&f=false
	,data.frame("Caranx hippos",'Crevalle jack',10)
	#https://en.wikipedia.org/wiki/King_mackerel
	,data.frame("Scomberomorus cavalla",'King mackerel',15)
	#http://myfwc.com/media/195458/spanish-mackerel.pdf
	#20^3.1373*.000193
	,data.frame("Scomberomorus maculatus",'Spanish mackerel',2.5)
	#http://www.tandfonline.com/doi/abs/10.1577/1548-8659%281962%2991%5B89%3ALWASLT%5D2.0.CO%3B2?journalCode=utaf20
	,data.frame("Trinectes maculatus",'Hogchoker',110/1000)
	#http://docserver.ingentaconnect.com/deliver/connect/umrsmas/00074977/v75n1/s6.pdf?expires=1468523196&id=88166574&titleid=10983&accname=University+of+Pennsylvania+Library&checksum=7CB64D8F5ACC7842E32F1ED74885792D
	#3.47*1e-6*500^3.21
	,data.frame("Paralichthys lethostigma",'Southern flounder',1.5)
	#https://www.flmnh.ufl.edu/fish/discover/species-profiles/carcharhinus-plumbeus
	,data.frame("Carcharhinus plumbeus",'Sandbar shark',65)
	#https://www.flmnh.ufl.edu/fish/discover/species-profiles/carcharhinus-brevipinna
	,data.frame("Carcharhinus brevipinna",'Spinner shark',56)
	#http://marinebio.org/species.asp?id=372
	,data.frame("Rhizoprionodon terraenovae",'Sharpnose shark',7.25)
)
weights<-do.call(rbind,lapply(tmp,function(x){colnames(x)<-c('species','name','weight');x}))
rownames(weights)<-weights$species


#downloaded files from SRA


#sample info from SRA
ncbi<-readLines('data/fish/biosample_result.txt')
sraLines<-grep('SRA:',ncbi)
idLines<-grep('original_id=',ncbi)
hostLines<-grep('/host=',ncbi)
lineDiff<-hostLines-idLines
if(any(lineDiff!=lineDiff[1]))stop(simpleError('Host and id lines not consistent'))

fishFirstLines<-sapply(list.files('data/fish','fastq$',full.names=TRUE),readLines,n=1)
#fishFirstLines<-sub('^@([^.]+)\\.[0-9].*$','\\1',fishFirstLines)
fishFirstLines<-sub('^[^ ]+ ([^_]+)_.*$','\\1',fishFirstLines)


info<-data.frame(
	'srs'=sub('.*SRA: ','',ncbi[sraLines])
	,'srr'=sub('.*/([^.]+).fastq','\\1',names(fishFirstLines))
	,'species'=sub('.*\\host="([^"]+)"','\\1',ncbi[hostLines])
	,'id'=sub('.*\\original_id="([^"]+)"','\\1',ncbi[idLines])
	,stringsAsFactors=FALSE
)

info$name<-weights[info$species,'name']
info$weight<-weights[info$species,'weight']
write.csv(info,'data/fish/sampleInfo.csv',row.names=FALSE)
info$file<-sprintf('data/fish/%s.fastq',info$srr)

dir.create('work/data/fish',showWarnings=FALSE)
info$common<-info$name
info$name<-info$file
write.csv(info[,c('name','species','common','weight')],'work/data/fish/info.csv')

allSeqs<-mclapply(info$file,function(x){
	tmp<-read.fastq(x,convert=TRUE)
	seqs<-filterReads(tmp$seq,tmp$qual,minLength=300,maxLength=375,minQual=10,maxBadQual=3)
	return(seqs)
},mc.cores=16)

runOtuForming(unlist(allSeqs),rep(info$file,sapply(allSeqs,length)),'work/data/fish')
