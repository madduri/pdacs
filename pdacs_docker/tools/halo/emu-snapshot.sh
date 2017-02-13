#!/bin/bash -l

if [ $# -ne 2 ]
then
  echo "Wrong number of arguments in $0: got $#, require 2"
  exit 1
fi

ofile=$2

omega_m=$(sqlite3 $1 "select value from metadata where name='Omega_m';")

hubble=$(sqlite3 $1 "select value from metadata where name='hubble';")

omega_b=$(sqlite3 $1 "select value from metadata where name='Omega_bar';")

n_s=$(sqlite3 $1 "select value from metadata where name='n_s';")

w=$(sqlite3 $1 "select value from metadata where name='w_de';")

sigma_8=$(sqlite3 $1 "select value from metadata where name='Sigma_8';")

a=$(sqlite3 $1 "select value from metadata where name='a_val';")
echo $a
z=`echo "scale=6;(1/$a)-(1)" | bc`
echo $z

tmp=`python -c "print float('$omega_b')"`
newomega=`echo "scale=6;(${hubble})*(${hubble})*(${omega_m})" | bc`
echo $newomega

tmp=`python -c "print float('$omega_b')"`

echo $tmp
newomegab=`echo "scale=6;(${hubble})*(${hubble})*(${tmp})" | bc`
echo $newomegab
echo $omega_b

echo $a $newomega $newomegab $n_s $sigma_8 $w $z 2
/global/project/projectdirs/hacc/PDACS/cosmic_emu/emu.exe emu.initial.output $newomega $newomegab $n_s $sigma_8 $w $z 2

status=$?
if [ $status -ne 0 ]
then
    echo emu failed with status $status
    exit 3
fi

#tr ' ' '\t' <initial.output >$ofile
( echo "#k P(k)" ; cat emu.initial.output ) > initial.output && rm emu.initial.output
sed '/^#/ d' initial.output | tr ' ' '\t' >$ofile



