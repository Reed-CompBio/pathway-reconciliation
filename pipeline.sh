#!/bin/bash

#this is just the list of commands needed to run the whole experiment. 


#first we initialize the graphlets directory if it doesn't exist

if [ ! -d graphlets ]
then
	mkdir graphlets
fi


#and even if it does, we populate it with any new dbs that aren't in it yet

for x in $(find networks/dbs -type d | cut -d "/" -f 3);
do
	if [ ! -d graphlets/$x ]
	then
		mkdir graphlets/"$x"
	fi
done

#now we populate those directories with graphlets and ghust coefficients, unless told not to

if [ "$1" = "generate" ]

then
	cd src/graphlets;

	for x in $(find ../../networks/dbs -type d | cut -d "/" -f 5);
	do
		python3 graphlet_count.py ../../networks/dbs/"$x"/ ../../graphlets/"$x"/
		python3 rho_coeff.py ../../graphlets/"$x"/ ../../graphlets/"$x"/
	done

	cd ../../
fi

#we also want to initialize an output directory to put our plots if it doesn't exist
if [ ! -d out ]
then
	mkdir out
fi


#now we can cluster dbs against one another in a pairwise fashion for instance. 

#for x in $(ls graphlets)
#do
#	for y in $(ls graphlets)
#	do
#		if [ ! "$x" = "$y" ]
#		then
#			cd src/clustering
#			python3 plot_dendrogram_colors.py ../../graphlets/"$x"/*.gcount ../../graphlets/"$y"/*.gcount ../../out/"$x"-vs-"$y"-dendrogram.pdf
#			cd ../../
#		fi
#	done
#done

#lets make one big dendrogram with every database
cd src/clustering
python3 plot_dendrogram_colors.py ../../graphlets/*/*.gcount ../../out/all-dbs-dendrogram.pdf
cd ../../

