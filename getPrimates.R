

weights<-c(
  "mantled howling monkey"=6, #https://en.wikipedia.org/wiki/Mantled_howler
  "red-shanked douc"=10, #https://en.wikipedia.org/wiki/Red-shanked_douc
  "Orangutan"=55, #https://en.wikipedia.org/wiki/Orangutan
  "Emperor tamarin"=.5, #https://en.wikipedia.org/wiki/Emperor_tamarin
  "Black-handed spider monkey"=7.5, #https://en.wikipedia.org/wiki/Geoffroy's_spider_monkey
  "White-faced saki"=1.8, #http://animaldiversity.org/accounts/Pithecia_pithecia/
  "Geoffroys tamarin"=.5, #https://en.wikipedia.org/wiki/Geoffroy's_tamarin
  "Blue-eyed black lemur"=1.8, #https://en.wikipedia.org/wiki/Blue-eyed_black_lemur
  "De Brazzas monkey"=5.5, #https://en.wikipedia.org/wiki/De_Brazza's_monkey
  "Western lowland gorilla"=125, #https://en.wikipedia.org/wiki/Western_lowland_gorilla
  "Sumatran Orangutan=77" #https://en.wikipedia.org/wiki/Sumatran_orangutan
)
info<-read.table('data/primates/map_final_analysis.txt',sep='\t',stringsAsFactors=FALSE,header=TRUE)
info$weight<-weights[info$SpeciesCommonNameSpecific]
info$species<-info$SpeciesScientificName
info$common<-info$SpeciesCommonNameSpecific
info$name<-rownames(info)

dir.create('work/data/primates',showWarnings=FALSE)
write.csv(info[,c('name','species','common','weight')],'work/data/primates/info.csv',row.names=FALSE)

seqs<-read.fa('data/primates/merged.fna.gz')
seqs$sample<-sub('_[0-9]+.*$','',seqs$name)
if(any(!seqs$sample %in% rownames(info)))stop(simpleError('Unknown sample found'))
seqs<-seqs[!grepl('[^ATCG]',seqs$seq),]

runOtuForming(substring(seqs$seq,1,min(nchar(seqs$seq))),seqs$sample,'work/data/primates')
