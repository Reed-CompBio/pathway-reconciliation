#!/bin/bash

upto=$1
path=${2?Error: No path given}
orca_path=$(pwd)/ORCA/orca/

cd $path

for i in *.orca
do	
	echo "running: ${orca_path}'orca.exe' node ${upto} ${i} ${i}'.ocount'";
	${orca_path}'orca.exe' node ${upto} ${i} ${i}'.ocount'
done
