#!/bin/bash -l

if [ $# -ne 5 ]
then
  echo "Wrong number of arguments in $0: got $#, require 5"
  exit 1
fi

echo $1 $2 $3 $4

ftype=$1
massBins=$2
metafile=$3
outputfile=$4
galaxy_dir=$5
tools_dir=`dirname ${galaxy_dir}`
#echo $tools_root


module unload pgi
module load gcc

foffile=$(sqlite3 $3 "select value from metadata where name='FOFPropertiesFile';")
sodfile=$(sqlite3 $3 "select value from metadata where name='SODPropertiesFile';")
cmfile=$(sqlite3 $3 "select value from metadata where name='SODPropertyBinsFile';")

#echo $foffile
#echo $sodfile
#echo $cmfile

ctmpfile=`echo "${cmfile}" | rev | cut -d/ -f1 |rev`
#cfile="${tools_dir}/working/${ctmpfile}"
cfile="${galaxy_dir}/database/tmp/${ctmpfile}"

if [ $ftype == 1 ]
then 
   tmpfile=`echo "${foffile}" | rev | cut -d/ -f1 |rev`
   value=${foffile}
elif [ $ftype == 2 ] || [ $ftype == 3 ]
then 
   tmpfile=`echo "${sodfile}" | rev | cut -d/ -f1 |rev`
   value=${sodfile}
fi

if [ $ftype == 3 ]
then
  #python /global/project/projectdirs/hacc/PDACS/pdacs-test/tools/halo/GenericIOPrint.py --input $cmfile --outFile $cfile --noRank 0
  echo "${tools_dir}/Frontend_reader/GenericIOPrint $cmfile > $cfile"
  ${tools_dir}/Frontend_reader/GenericIOPrint $cmfile > $cfile
fi

#propfile="${tools_dir}/working/${tmpfile}"
propfile="${galaxy_dir}/database/tmp/${tmpfile}"

#python /global/project/projectdirs/hacc/PDACS/pdacs-test/tools/halo/GenericIOPrint.py --input $value --outFile $propfile --noRank 0
${tools_dir}/Frontend_reader/GenericIOPrint $value > $propfile

#echo `pwd`
python ${galaxy_dir}/tools/halo/driver.py -j $ftype -m $massBins -i $metafile -p $propfile -c $cfile -o $outputfile

#rm "${propfile}"

