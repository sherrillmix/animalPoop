#!/bin/bash
##File name: swarmData.bash
##Creation date: Aug 24, 2015
##Last modified: Thu May 26, 2016  04:00PM
##Created by: scott
##Summary: Swarm to make OTUs

cd data/mgrast

echo Starting filtering fas
for ii in *.fa.gz;do 
	echo $ii
	sem -j8  Rscript ../../uniquify.R $ii ${ii%.fa.gz}_uniq.fa 250 350
done
sem --wait
echo Done filtering fas

echo Starting filtering fastqs
for ii in *.fastq.gz;do 
	echo $ii
	sem -j8  Rscript ../../uniquify.R $ii ${ii%.fastq.gz}_uniq.fa 250 350
done
sem --wait
echo Done filtering fastqs

echo Starting swarm
for ii in *_uniq.fa;do 
	echo $ii
	sem -j8  ~/installs/swarm/swarm $ii -o ${ii%_uniq.fa}.out
done
sem --wait
echo Done swarming

echo All done
