#!/bin/bash -l

if [ $# -ne 19 ]
then
  echo "Wrong number of arguments in $0: got $#, require 19"
  exit 1
fi

ftype=$1
mBins=$2
nFiles=$3
zSeek=$4
mNum=$5
hubble=$6
omegaM=$7
ns=$8
sigma8=$9
w0=${10}
boxsize=${11}
particlenum=${12}
minparticle=${13}
delta=${14}
haloFile=${15}
inputFile=${16}
tkFile=${17}
dirPath=${18}
outfile=${19}

# Run the actual program.
if [ -f $outfile ]
then
  rm $outfile
fi

echo "$outfile"
python /global/project/projectdirs/hacc/PDACS/galaxy_dev_ravi/tools/halo/MassFunction.py -j $ftype -m $mBins -n $nFiles -z $zSeek -l $mNum -e $hubble -g $omegaM -a $ns -s $sigma8 -w $w0 -b $boxsize -p $particlenum -r $minparticle -d $delta -u $haloFile -i $inputFile -t $tkFile -c $dirPath -o $outfile 

status=$?
if [ $status -ne 0 ]
then
    echo MassFunction failed with status $status
    exit 3
fi

#This needs to be executed for all
sqlite3 $outfile <<EOF
create table metadata (
  name TEXT,
  value NUMERIC);

insert into metadata ( name, value ) values ( "program", "MassFunction" );
insert into metadata ( name, value ) values ( "mBins", "${mBins}" );
insert into metadata ( name, value ) values ( "nFiles", "${nFiles}" );
insert into metadata ( name, value ) values ( "zSeek", "${zSeek}" );
insert into metadata ( name, value ) values ( "mNum", "${mNum}" );
insert into metadata ( name, value ) values ( "hubble", "${hubble}" );
insert into metadata ( name, value ) values ( "omegaM", "${omegaM}" );
insert into metadata ( name, value ) values ( "ns", "${ns}" );
insert into metadata ( name, value ) values ( "sigma8", "${sigma8}" );
insert into metadata ( name, value ) values ( "w0", "${w0}" );
insert into metadata ( name, value ) values ( "boxsize", "${boxsize}" );
insert into metadata ( name, value ) values ( "particlenum", "${particlenum}" );
insert into metadata ( name, value ) values ( "minparticle", "${minparticle}" );
insert into metadata ( name, value ) values ( "delta", "${delta}" );
insert into metadata ( name, value ) values ( "haloFile", "${haloFile}" );
insert into metadata ( name, value ) values ( "inputFile", "${inputFile}" );
insert into metadata ( name, value ) values ( "tkFile", "${tkFile}" );
insert into metadata ( name, value ) values ( "dirPath", "${dirPath}" );
insert into metadata ( name, value ) values ( "outfile", "${outfile}" );
EOF

