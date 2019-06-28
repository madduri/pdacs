#!/bin/bash

echo $1 $2 $3 $4 $5 $6 | tr ' ' \\n > cm.local.param.ini
# Run the actual program.
/project/projectdirs/hacc/PDACS/cm/cM cm.local.param.ini cm.initial.output
#grep -v '^#' initial.output > stripped.output
( echo "logmass	conc" ; cat cm.initial.output ) > initial.output && rm cm.initial.output
sed '/^#/ d' initial.output | tr ' ' '\t' >$7
rm initial.output
