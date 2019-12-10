#!/bin/bash

# need a better way, but right now assumption is that files reside under this dirtree /project/projectdirs/hacc/PDACS/Coyote
#dirname=`find /project/projectdirs/hacc/PDACS -name "Coyote" 2>/dev/null -print`

#change to the simulation directory
#echo ${dirname}
#cd ${dirname}
cd  /project/projectdirs/hacc/PDACS/Coyote
# find all snapshot files in the simulation directory 
find . -type f -name "snapshot_*.*.0" 2>/dev/null > ../Catalog/tmp.dat

cd ../Catalog
export IFS="/"
printf "Simulation\tType\tModel\tSize\tRealization\tSnapshot\tFull-Path\n" > coyote.txt
cat tmp.dat | while read a b c d e f; do printf "Coyote\t$b\t$c\t$d\t$e\t${f%??}\t$c/$d/$e/${f%??}\n" >> coyote.txt; done

