#!/bin/bash
##File name: swarmData.bash
##Creation date: Aug 24, 2015
##Last modified: Tue Aug 25, 2015  10:00AM
##Created by: scott
##Summary: Swarm to make OTUs

cd data/bushman

if [ ! -f uniq.fa.gz ];then
	echo "Finding unique reads in Bushman data"
	Rscript ../../uniquifyBushman.R
	echo "Done finding unique reads in Bushman data"
fi

gunzip uniq.fa.gz

echo Starting swarm
~/installs/swarm/swarm uniq.fa -o uniq.out
echo Done swarming

gzip uniq.fa

echo All done
