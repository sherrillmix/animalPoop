#!/bin/bash

cd data/bushman/split/

echo Starting filtering
for ii in *.fa.gz;do 
	echo $ii
	sem -j8 Rscript ../../../uniquify.R $ii ${ii%.fa.gz}_uniq.fa 250 350 #sem from apt-get install parallel 
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
