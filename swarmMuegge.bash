#!/bin/bash
##File name: swarmData.bash
##Creation date: Aug 24, 2015
##Last modified: Tue May 31, 2016  09:00AM
##Created by: scott
##Summary: Swarm to make OTUs

cd data/muegge

echo Starting filtering
for ii in *.fastq.gz;do 
	echo $ii
	sem -j8  Rscript ../uniquify.R $ii ${ii%.fastq.gz}.fa #sem from apt-get install parallel 
done
sem --wait
echo Done filtering

echo Starting swarm
for ii in *.fa;do 
	echo $ii
	sem -j8  ~/installs/swarm/swarm $ii -o ${ii%.fa}.out
done
sem --wait
echo Done swarming

echo All done
