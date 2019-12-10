#!/bin/bash

# need a better way, but right now assumption is that files reside under this dirtree /project/projectdirs/hacc/PDACS/Coyote
dirname=`find /project/projectdirs/hacc/PDACS/Coyote/Grid/$1/$2 -name "$3" 2>/dev/null -print`

#change to the simulation directory
echo ${dirname}
cd ${dirname}

# find all snapshot files in the simulation directory 
find . -type f -name "snapshot_*.0" 2>/dev/null > ../../../../../Catalog/tmp1.dat

export IFS="="
while read key val;
do
if [[ $key == "box_size" ]]
then
box_size=$val;
fi
if [[ $key == "seed" ]]
then
seed=$val;
fi
if [[ $key == "z_in" ]]
then
z_in=$val;
fi
done < params.ini
echo ${box_size}
echo ${seed}
echo ${z_in}

while read val;
do
val=`echo "$val"|cut -d "/" -f2`
sqlite3 /project/projectdirs/hacc/PDACS/Catalog/$1 <<EOF
insert into SnapshotObj ( Snapshot, SnapshotPath, BoxSizeName, BoxSize, Realization, Seed, zin ) values ( "${val//.0/}", "/project/projectdirs/hacc/PDACS/Coyote/Grid/$1/$2/$3/${val//.0/}", "$2", "${box_size}", "$3", "${seed}", "${z_in}");
EOF
done < ../../../../../Catalog/tmp1.dat

