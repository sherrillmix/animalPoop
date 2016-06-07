#!/bin/bash

cd data/anteater/

echo Starting filtering
for ii in *.fastq.gz;do 
	echo $ii
	sem -j8 Rscript ../../uniquify.R $ii ${ii%.fastq.gz}_uniq.fa #sem from apt-get install parallel 
done
sem --wait
echo Done filtering

echo Starting swarm
for ii in *_uniq.fa;do 
	echo $ii
	sem -j8 ~/installs/swarm/swarm $ii -o ${ii%.fa}.out
done
sem --wait
echo Done swarming

echo All done
