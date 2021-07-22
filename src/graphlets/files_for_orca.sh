#!/bin/bash

path=${1?Error: No path given}

cd $path

for i in *network.txt
do
	m=$(wc -l < $i)
	n=$(grep -Eo '[0-9]+' $i | sort -rn | head -n 1)
	n=$(($n + 1))
	echo $n $m > ${i}'.orca'
	cat $i >> ${i}'.orca'
done
