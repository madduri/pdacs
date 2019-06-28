#!/bin/bash

if [ $# -ne 8 ]
then
  echo "Wrong number of arguments in $0: got $#, require 8"
  exit 1
fi

ofile=$1
newomega=$2
newomegab=$3
n_s=$4
sigma_8=$6
w=$5
z=$7
otype=$8

touch initial.output
if [ -f initial.output ]
then
  echo -n "" > initial.output
fi

ls initial.output
/global/project/projectdirs/hacc/PDACS/cosmic_emu/emu.exe emu.initial.output $newomega $newomegab $n_s $sigma_8 $w $z 2
status=$?
if [ $status -ne 0 ]
then
    echo emu failed with status $status
    exit 3
fi

( echo "k Pk" ; cat emu.initial.output ) > initial.output && rm emu.initial.output
sed '/^#/ d' initial.output | tr ' ' '\t' >$ofile
rm initial.output

