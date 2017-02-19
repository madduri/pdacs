import sys
import os
from subprocess import call
import optparse
import string
from time import time, sleep
import json
import shutil
import sqlite3 as lite 
import fileinput

#mpirun ... xi.out ngrid inputfile nfiles outputfile boxsize [line of sight w omega_m]

#mpirun -n %(np)d %(toolsDir)s/xi.out %(ngrid)s %(inputfile)s %(nfiles)s %(outfile)s %(boxsize)s %(line)s %(w)s %(omega_m)s

def run():    
    parser = optparse.OptionParser()
    parser.add_option('--inputfile', action="store", dest="inputfile")
    parser.add_option('--ngrid', action="store", dest="ngrid")
    parser.add_option('--outfile_2d', action="store", dest="outfile_2d")
    parser.add_option('--outfile_multipole', action="store", dest="outfile_multipole")
    parser.add_option('--line', action="store", dest="line")
    parser.add_option('--numNodes', action="store", dest="numNodes")
    parser.add_option('--numProcs', action="store", dest="numProcs")    
    parser.add_option('--exeTime', action="store", dest="exeTime")    
    parser.add_option('--exeQueue', action="store", dest="exeQueue")    

    (options, args) = parser.parse_args()
    
    nodes = options.numNodes
    procs = options.numProcs

    toolsDir = "/global/project/projectdirs/hacc/PDACS/JK_Tools/"
    workingDir = "./"

    #get the name of the file - used to make the temp output files.
    workingName = options.outfile_2d.split("/")
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
    cur.execute("select value from metadata where name='Omega_m'")
    row = cur.fetchone()
    omega_m = float(row[0])
    cur.execute("select value from metadata where name='w_de'")
    row = cur.fetchone()
    w = float(row[0])
    cur.execute("select value from metadata where name='Snapshot'")
    con.commit()
    row = cur.fetchone()
    snapshotname = row[0]

    infile = "/global/project/projectdirs/hacc/PDACS/Coyote/Grid/" + snapshotname
    #print infile
    
    #write the pbs file to execute on Carver
    #pbsPath = create_pbs_file(nodes, procs, options.exeQueue, options.exeTime, toolsDir, workingDir, workingName, options.ngrid, infile, nfiles, outFile, boxsize, options.line, w, omega_m)
    pbsCmd = "srun -n %s %s/xi.out %s %s %s %s %s %s %s %s" % (nodes, toolsDir, options.ngrid, infile, nfiles, outFile, boxsize, options.line, w, omega_m)
    
    os.system(pbsCmd)
    
    #TODO what is this file extension when made for these two?

    cmd = "sed -e 's/^[ \t]*//;s/[ \t]*$//' %s | tr -s ' ' '\t' > %s" %(outFile + ".2D", options.outfile_2d)
    os.system(cmd)
    headers = '#r_perp r_parallel mu xi(r_perp, r_para)'.split()
    for line in fileinput.input([options.outfile_2d], inplace=True):
      if fileinput.isfirstline():
        print '\t'.join(headers)
      if len(line.strip())>0:
        print line

    cmd = "sed -e 's/^[ \t]*//;s/[ \t]*$//' %s | tr -s ' ' '\t' > %s" %(outFile + ".multipole", options.outfile_multipole)
    os.system(cmd)
    headers = '#r monopole quadrupole'.split()
    for line in fileinput.input([options.outfile_2d], inplace=True):
      if fileinput.isfirstline():
        print '\t'.join(headers)
      if len(line.strip())>0:
        print line

if __name__ == '__main__':
    run()
