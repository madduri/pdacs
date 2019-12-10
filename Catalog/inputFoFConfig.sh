#!/bin/bash

# need a better way, but right now assumption is that files reside under this dirtree /project/projectdirs/hacc/PDACS/Coyote
dirname=`find /project/projectdirs/hacc/PDACS -name "Coyote" 2>/dev/null -print`

#change to the simulation directory
echo ${dirname}
cd ${dirname}

# find all snapshot files in the simulation directory 
find . -type f -name "input.M001" 2>/dev/null > ../Catalog/tmp1.dat

cd ../Catalog
#export IFS="/"
#printf "Simulation\tType\tModel\tSize\tRealization\tSnapshot\tFull-Path\n" > coyote.txt
#cat tmp.dat | while read a b c d e f; do printf "Coyote\t$b\t$c\t$d\t$e\t${f//.0/}\t$c/$d/$e/${f//.0/}\n" >> coyote.txt; done

