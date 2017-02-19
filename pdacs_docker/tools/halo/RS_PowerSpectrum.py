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
#mpirun ... powerspec.out ngrid inputfile nfiles outputfile boxsize [line of sight w omega_m]


def run():    
    parser = optparse.OptionParser()
    parser.add_option('--inputfile', action="store", dest="inputfile")
    parser.add_option('--ngrid', action="store", dest="ngrid")
    parser.add_option('--exeTime', action="store", dest="exeTime")
    parser.add_option('--exeQueue', action="store", dest="exeQueue")
    parser.add_option('--outfile_2d', action="store", dest="outfile_2d")
    parser.add_option('--outfile_multipole', action="store", dest="outfile_multipole")
#    parser.add_option('--boxsize', action="store", dest="boxsize")
    parser.add_option('--line', action="store", dest="line")
#    parser.add_option('--w', action="store", dest="w")
#    parser.add_option('--omega_m', action="store", dest="omega_m")
    parser.add_option('--numNodes', action="store", dest="numNodes")
    parser.add_option('--numProcs', action="store", dest="numProcs")    

    (options, args) = parser.parse_args()
    
    nodes = options.numNodes
    procs = options.numProcs

    numnodes = int(nodes)
    numprocs = int(procs)
    np = numnodes*numprocs

    toolsDir = "/global/project/projectdirs/hacc/PDACS/JK_Tools/"
    workingDir = "./"

   # inputArg = options.inputfile
   # split = os.path.splitext(options.inputfile)

    #get the name of the file - used to make the temp output files.
    workingName = options.outfile_2d.split("/")
    workingName = workingName[len(workingName)-1]
    #outFile = "rs_pow.out"
    #workingDir + workingName
    #set the path to create files
    #if not os.path.exists(workingDir):
    #    os.makedirs(workingDir)

    outFile = workingDir + workingName + ".out"
    #write the pbs file to execute on Carver
    #os.chmod(outFile,S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH|S_IWOTH)

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

#    print options.outfile
    infile = "/global/project/projectdirs/hacc/PDACS/Coyote/Grid/" + snapshotname
#    print infile
 
    #pbsPath = create_pbs_file(nodes, procs, options.exeQueue, options.exeTime, toolsDir, workingDir, workingName, options.ngrid, infile, nfiles, outFile, boxsize, options.line, w, omega_m)
    pbsCmd = "module unload pgi; module load gcc; module load fftw; module load gsl; module load sqlite; srun -n %s %s/powerspec.out %s %s %s %s %s %s %s %s" % (numnodes, toolsDir, options.ngrid, infile, nfiles, outFile, boxsize, options.line, w, omega_m)
    print "CMD: " + pbsCmd
    print "CWD: " + os.getcwd()
    os.system(pbsCmd)

    #shutil.copyfile(outFile, options.outputfile)
    cmd = "sed -e 's/^[ \t]*//;s/[ \t]*$//' %s | tr -s ' ' '\t' > %s" %(outFile+".2D", options.outfile_2d)
    os.system(cmd)
    headers = '#k mu P(k,mu)'.split()
    for line in fileinput.input([options.outfile_2d], inplace=True):
      if fileinput.isfirstline():
        print '\t'.join(headers)
      if len(line.strip())>0:
        print line
    cmd = "sed -e 's/^[ \t]*//;s/[ \t]*$//' %s | tr -s ' ' '\t' > %s" %(outFile+".multipole", options.outfile_multipole)
    os.system(cmd)
    headers = '#k monopole quadrupole'.split()
    for line in fileinput.input([options.outfile_multipole], inplace=True):
      if fileinput.isfirstline():
        print '\t'.join(headers)
      if len(line.strip())>0:
        print line
    #os.remove(outFile+".2D")
    os.remove(outFile+".multipole")

if __name__ == '__main__':
    run()
