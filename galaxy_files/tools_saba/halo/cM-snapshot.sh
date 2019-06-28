#!/bin/bash

omega_m=$(sqlite3 $1 "select value from metadata where name='Omega_m';")

omega_b=$(sqlite3 $1 "select value from metadata where name='Omega_bar';")

n_s=$(sqlite3 $1 "select value from metadata where name='n_s';")

w=$(sqlite3 $1 "select value from metadata where name='w_de';")

sigma_8=$(sqlite3 $1 "select value from metadata where name='Sigma_8';")

a=$(sqlite3 $1 "select value from metadata where name='a_val';")

hubble=$(sqlite3 $1 "select value from metadata where name='hubble';")

z=`echo "(1/$a)-(1)" | bc`
#echo "red shift is $z"

hsqr=`echo "($hubble)*($hubble)" | bc`
#echo "hsqr: $hsqr"
newomega=`echo "($hsqr) * ($omega_m)" | bc`
#echo "Omega_m: $newomega"

tmp=`python -c "print float('$omega_b')"`
newomegab=`echo "($hubble)*($hubble)*($tmp)" | bc`

echo $newomega $newomegab $n_s $w $sigma_8 $z | tr ' ' \\n > local.param.ini

#echo initial.output
#echo $6
# Run the actual program.
/project/projectdirs/hacc/PDACS/cm/cM local.param.ini cm.initial.output
#grep -v '^#' initial.output > stripped.output
( echo "logmass conc" ; cat cm.initial.output ) > initial.output && rm cm.initial.output
sed '/^#/ d' initial.output | tr ' ' '\t' >$2
rm initial.output
