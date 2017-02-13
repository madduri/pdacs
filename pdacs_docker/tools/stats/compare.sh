#!/bin/bash

#if [ $# -ne 2 ]
#then
#  echo "Wrong number of arguments in $0: got $#, require 2"
#  exit 1
#fi

echo $@
numargs=$#
dbname=$2
#numcols=`awk -F"\t" '{print NF}' $infile`

sqlite3 $dbname "create table metadata(tabname TEXT);"

for (( c=1; c<=${numargs}; c=c+3 ))
do
   filename=$1
   shift
   shift
   shift
   tablename=`echo "${filename}" | rev | cut -d/ -f1 |rev | cut -d. -f1`
   echo $tablename
#tablename="tst"
#dbname=$2 
#echo $dbname
   sed 's/\t/,/g' $filename > output_file
   list=$(cat output_file | sed -n '1p' | sed 's/,/ NUMERIC, /g')
   echo $list
   cat output_file | sed '1d' | sed 's/,/|/g' | sed '/^\s*$/d' > temp
   sqlite3 $dbname "create table $tablename($list NUMERIC);"
   sqlite3 $dbname ".import temp $tablename"
   rm temp
done
