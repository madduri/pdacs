import sys, pwd, os, optparse, fileinput
from subprocess import call
from time import time, sleep
import sqlite3 as lite

#mpirun ... powerspec.out ngrid inputfile nfiles outputfile boxsize [line of sight w omega_m]

def run():    
    parser = optparse.OptionParser()
    parser.add_option('--inputfile', action="store", dest="inputfile")
    parser.add_option('--ngrid', action="store", dest="ngrid")
    parser.add_option('--outfile', action="store", dest="outfile")
    parser.add_option('--numNodes', action="store", dest="numNodes")
    parser.add_option('--numProcs', action="store", dest="numProcs")    
    parser.add_option('--exeTime', action="store", dest="exeTime")    
    parser.add_option('--exeQueue', action="store", dest="exeQueue")    

    (options, args) = parser.parse_args()
    
    nodes = options.numNodes
    procs = options.numProcs

    ngrid = options.ngrid
    numnodes = int(nodes)
    numprocs = int(procs)

    toolsDir = "/global/project/projectdirs/hacc/PDACS/JK_Tools"

    #get the name of the file - used to make the temp output files.
    workingName = options.outfile.split("/")
    workingName = workingName[len(workingName)-1]

    outFile = "./" + workingName + ".out"
 
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
    
    #print options.outfile
    inFile = "/global/project/projectdirs/hacc/PDACS/Coyote/Grid/" + snapshotname
    # print inFile

    #write the pbs file to execute on Carver
    #pbsCmd = "module unload pgi; module load gcc; module load fftw; module load gsl; module load sqlite; mpirun -n %d %s/powerspec.out %s %s %f %s %f" % (numnodes*numprocs, toolsDir, ngrid, options.inputfile, nfiles, outFile, boxsize)
    pbsCmd = "module load fftw; module load gsl; module load sqlite; mpirun -n %d %s/powerspec.out %s %s %f %s %f" % (numnodes*numprocs, toolsDir, ngrid, inFile, nfiles, outFile, boxsize)
    os.system(pbsCmd)
    print "CMD: " + pbsCmd
    print "CWD: " + os.getcwd()

    cmd = "sed 's/^[ \t]*//' %s | tr -s ' ' '\t' > %s" %(outFile, options.outfile)
    os.system(cmd)
    headers = '#K PK'.split()
    for line in fileinput.input([options.outfile], inplace=True):
      if fileinput.isfirstline():
        print '\t'.join(headers)
      if len(line.strip())>0:
        print line

    os.remove(outFile)
if __name__ == '__main__':
    run()
