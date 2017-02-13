#!/usr/bin/python

import sys, os
from subprocess import call
import optparse
import string
from time import time, sleep
import json
import shutil
import sqlite3 as lite
import fileinput
#mpirun ... xi.out ngrid inputfile nfiles outputfile boxsize [line of sight w omega_m]

def run():    
    parser = optparse.OptionParser()
    parser.add_option('--inputfile', action="store", dest="inputfile")
    parser.add_option('--ngrid', action="store", dest="ngrid")
    parser.add_option('--exeTime', action="store", dest="exeTime")
    parser.add_option('--exeQueue', action="store", dest="exeQueue")
    parser.add_option('--outfile', action="store", dest="outfile")
    parser.add_option('--outfile_pk', action="store", dest="outfile_pk")
#    parser.add_option('--boxsize', action="store", dest="boxsize")
    parser.add_option('--numNodes', action="store", dest="numNodes")
    parser.add_option('--numProcs', action="store", dest="numProcs")    

    (options, args) = parser.parse_args()
    
    nodes = options.numNodes
    procs = options.numProcs

    toolsDir = "/global/project/projectdirs/hacc/PDACS/JK_Tools/"
    #workingDir = "/global/project/projectdirs/hacc/PDACS/working/"
    workingDir = "./"

   # inputArg = options.inputfile
   # split = os.path.splitext(options.inputfile)

    #get the name of the file - used to make the temp output files.
    workingName = options.outfile.split("/")
    workingName = workingName[len(workingName)-1]
    
    outFile = workingDir + workingName + ".out"

    con = lite.connect(options.inputfile)
    cur = con.cursor()
    cur.execute("select value from metadata where name='box_size [Mpc/h]'")
    row = cur.fetchone()
    boxsize = float(row[0])
    cur.execute("select value from metadata where name='numFiles'")
    row = cur.fetchone()
    nfiles = float(row[0])
    cur.execute("select value from metadata where name='Snapshot'")
    con.commit()
    row = cur.fetchone()
    snapshotname = row[0]
    
    infile = "/global/project/projectdirs/hacc/PDACS/Coyote/Grid/" + snapshotname 
    #write the pbs file to execute on Carver
    #pbsPath = create_pbs_file(nodes, procs, options.exeQueue, options.exeTime, toolsDir, workingDir, workingName, options.ngrid, infile, nfiles, outFile, boxsize)
    pbsCmd = "srun -n %s %s/xi.out %s %s %s %s %s" % (nodes, toolsDir, options.ngrid, infile, nfiles, outFile, boxsize)
    
    os.system(pbsCmd)
    
    #preprocess the output files
    cmd = "sed -e 's/^[ \t]*//;s/[ \t]*$//' %s | tr -s ' ' '\t' > %s" %(outFile, options.outfile)
    os.system(cmd)
    headers = '#k mu'.split()
    for line in fileinput.input([options.outfile], inplace=True):
      if fileinput.isfirstline():
        print '\t'.join(headers)
      if len(line.strip())>0:
        print line

    cmd = "sed -e 's/^[ \t]*//;s/[ \t]*$//' %s | tr -s ' ' '\t' > %s" %(outFile + ".pk", options.outfile_pk)
    os.system(cmd)
    headers = '#k pk'.split()
    for line in fileinput.input([options.outfile_pk], inplace=True):
      if fileinput.isfirstline():
        print '\t'.join(headers)
      if len(line.strip())>0:
        print line


if __name__ == '__main__':
    run()
