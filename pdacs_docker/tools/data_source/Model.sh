#!/bin/bash

# need a better way, but right now assumption is that files reside under this dirtree /global/project/projectdirs/hacc/PDACS/Coyote
#dirname=`find /global/project/projectdirs/hacc/PDACS -name "Coyote" 2>/dev/null -print`

#change to the simulation directory
#echo ${dirname}
#cd ${dirname}

# find all snapshot files in the simulation directory 
#find . -type f -name "params.ini" 2>/dev/null > ../Catalog/params-loc.dat

#export IFS="/"

#printf "ModelGrid\tType\tModel\tSize\tRealization\n" > coyote.txt
#cat ../Catalog/params-loc.dat | while read a b c d e; do printf "Coyote.$a.$b.$c.$d.$e\n" >> ../Catalog/params-name.dat; done

#while read $simulation $model $grid $L $G; 
#do 
        #name="$simulation.$model.$grid.$L.$G.$param\n";
        #path="$simulation/$model/$grid/$L/$G/$param\n";
	#echo -e "$name $path";
#        echo -e "$simulation $model $grid $L $G\n";
#	export IFS="="
        #while read key val;
	#do
	#	#echo -e "$key\t$val\n";
	#done < "$line1"
#done < ../Catalog/params-loc.dat

# Now create an SQLite3 database

sqlite3 $5 <<EOF
create table SnapshotObj (
 BoxSizeName TEXT,
 BoxSize REAL,
 Realization TEXT,
 Seed REAL,
 zin REAL,
 Snapshot TEXT,
 SnapshotPath TEXT);

create table metadata (
  name TEXT,
  value NUMERIC);

insert into metadata (name, value) values ( "Model", "$1");
insert into metadata (name, value) values ( "Size [Mpc]", "$2");
insert into metadata (name, value) values ( "Realization", "$3");
insert into metadata (name, value) values ( "Snapshot", "$4");

EOF

export IFS="="
while read key val;
do
#if [[ $key != *//* ]] && [[ $key != "" ]] && [[ $key != "box_size" ]] && [[ $key != "seed" ]] && [[ $key != "z_in" ]] && [[ $key != "PrintFormat" ]] 
if [[ $key != *//* ]] && [[ $key != "" ]]
then
   if [[ $key == "box_size" ]] 
     then 
	key="box_size [Mpc/h]"
   fi
   sqlite3 $5 <<EOF
   insert into metadata ( name, value ) values ( "$key", "$val" );
EOF
fi

done < /global/project/projectdirs/hacc/PDACS/Coyote/Grid/$1/$2/$3/params.ini 

numfiles=`ls -l /global/project/projectdirs/hacc/PDACS/Coyote/Grid/${4}* | wc -l`

#not a very good way but just trying to get line number
snapshot=`echo "$4" | cut -d "/" -f4`
a_val=`echo "$snapshot" | cut -d "_" -f2`
#linenum=`echo $linenum|sed 's/^0*//'`
#linenum=$((linenum + 1))
#echo ${linenum}
#a_val=`sed -n ${linenum}p /global/project/projectdirs/hacc/PDACS/Coyote/Grid/$1/$2/$3/aout.dat`

#if [[ ${a_val} == "" ]]
#then
#  linenum=$((linenum - 1))
#  a_val=`sed -n ${linenum}p /global/project/projectdirs/hacc/PDACS/Coyote/Grid/$1/$2/$3/aout.dat`
#fi

sqlite3 $5 <<EOF
insert into metadata ( name, value ) values ( "numFiles", "${numfiles}" );
insert into metadata ( name, value ) values ( "a_val", "${a_val}" );
EOF

#done < /global/project/projectdirs/hacc/PDACS/Coyote/Grid/$1/$2/$3/aout.dat

#export IFS="/"
#printf "Simulation\tType\tModel\tSize\tRealization\tSnapshot\tFull-Path\n" > coyote.txt
#cat tmp.dat | while read a b c d e f; do printf "Coyote\t$b\t$c\t$d\t$e\t${f//.0/}\t$c/$d/$e/${f//.0/}\n" >> coyote.txt; done

