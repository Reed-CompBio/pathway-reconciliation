#!/bin/bash

#this is just the list of commands needed to run the whole experiment. 


#first we initialize the graphlets directory if it doesn't exist

if [ ! -d graphlets ]
then
	mkdir graphlets
	mkdir graphlets/dbs
	mkdir graphlets/null-models
fi


#and even if it does, we populate it with any new dbs that aren't in it yet

for x in $(find networks/dbs -type d | cut -d "/" -f 3);
do
	if [ ! -d graphlets/dbs/$x ]
	then
		mkdir graphlets/dbs/"$x"
	fi

	if [ ! -d graphlets/null-models/$x ]
	then
		mkdir graphlets/null-models/"$x"
	fi
done

#now we populate the dbs directories with graphlets and ghust coefficients, unless told not to

if [ ! "$1" = "False" ]

then
	cd src/graphlets;

	for x in $(find ../../networks/dbs -type d | cut -d "/" -f 5);
	do
		python3 graphlet_count.py ../../networks/dbs/"$x"/ ../../graphlets/dbs/"$x"/
		python3 rho_coeff.py ../../graphlets/dbs/"$x"/ ../../graphlets/dbs/"$x"/
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
python3 plot_dendrogram_colors.py ../../graphlets/dbs/*/*.gcount ../../out/all-dbs-dendrogram.pdf
cd ../../


#next lets handle making random networks.

#we have two random network models. (1) random-empirical, (2) random-rewiring.


if [ ! -d networks/null-models/random-empirical ]
then
	mkdir networks/null-models/random-empirical

fi

if [ ! -d networks/null-models/random-rewiring ]
then
	mkdir networks/null-models/random-rewiring
fi

#populate both directories

for x in $(find networks/dbs -type d | cut -d "/" -f 3);
do
	if [ ! -d networks/null-models/random-empirical/$x ]
	then
		mkdir networks/null-models/random-empirical/"$x"
	fi

	if [ ! -d networks/null-models/random-rewiring/$x ]
	then
		mkdir networks/null-models/random-rewiring/"$x"
	fi
done


#if we want to generate null models, then we will do so.

if [ ! "$2" = "False" ]
then
	cd src/null-models
	for x in $(find ../../networks/dbs -type d | cut -d "/" -f 5);
	do
		echo $x
		python3 random-empirical.py ../../networks/interactomes/All_Pathway_Commons.txt ../../graphlets/dbs/"$x" 10 ../../networks/null-models/random-empirical/"$x";
	done
	cd ../../
fi


