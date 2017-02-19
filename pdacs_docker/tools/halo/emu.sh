#!/bin/bash -l

if [ $# -ne 8 ]
then
  echo "Wrong number of arguments in $0: got $#, require 8"
  exit 1
fi

ofile=$1
omega_m=$2
omega_b=$3
n_s=$4
sigma_8=$5
w=$6
z=$7
otype=$8

# Run the actual program.
if [ -f $ofile ]
then
  rm $ofile
fi

echo $ofile
echo $omaga_m
echo $omega_b
echo $n_s
echo $sigma_8

echo /global/project/projectdirs/hacc/PDACS/cosmic_emu/emu.exe initial.output $omega_m $omega_b $n_s $w $sigma_8 $z 2
/global/project/projectdirs/hacc/PDACS/cosmic_emu/emu.exe initial.output $omega_m $omega_b $n_s $w $sigma_8 $z 2

#tr ' ' '\t' <initial.output >$ofile

sed '/^#/ d' initial.output | tr ' ' '\t' >$ofile
sed "1i#K\tPk" $ofile

#grep -v '^#' initial.output > stripped.output

#if [ ! -f stripped.output ]
#then
#  echo Failed to create stripped.output
#  exit 2
#fi

# Now create an SQLite3 database
#sqlite3 $ofile <<EOF
#create table emu (
# k    NUMERIC,
# Pk   NUMERIC);

#create table metadata (
#  name TEXT,
#  value NUMERIC);

#.separator " "
#.import stripped.output emu
#insert into metadata ( name, value ) values ( "program", "emu" );
#insert into metadata ( name, value ) values ( "Omega_m h^2", "${omega_m}" );
#insert into metadata ( name, value ) values ( "Omega_b h^2", "${omega_b}" );
#insert into metadata ( name, value ) values ( "n_s", "${n_s}" );
#insert into metadata ( name, value ) values ( "w", "${w}" );
#insert into metadata ( name, value ) values ( "sigma_8", "${sigma_8}" );
#insert into metadata ( name, value ) values ( "z", "${z}" );
#EOF
