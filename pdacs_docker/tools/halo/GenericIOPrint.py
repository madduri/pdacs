import sys
import os
from subprocess import call
import optparse
import httplib2
import urllib, logging
import string
from time import time, sleep
import json
import shutil
import NewtHelper

def create_pbs_file(workingDir, flags, inputArg, outputFile, workingName):    
    pbsText = """
#PBS -q serial
#PBS -A hacc
##PBS -l nodes=1:ppn=1
#PBS -l walltime=00:5:00
#PBS -N ioprint
#PBS -e ioprint.$PBS_JOBID.err
#PBS -o ioprint.$PBS_JOBID.out
#PBS -V

module unload pgi
module load gcc

cd %s

/global/project/projectdirs/hacc/PDACS/HaloFind_v3/GenericIOPrint%s%s > %s
chmod 777 %s
chgrp hacc %s
""" % (workingDir, flags, inputArg, outputFile, outputFile, outputFile)
    pbsPath = workingDir + workingName + ".pbs"
    pbsFile = open(pbsPath, 'w')
    pbsFile.write(pbsText)
    pbsFile.close()
    print 'Created PBS file.'
    return pbsPath
    

def run():
    parser = optparse.OptionParser()
    parser.add_option('--input', action="store", dest="input")
    parser.add_option('--outFile', action="store", dest="outFile")
    parser.add_option('--noRank', action="store", dest="noRank")
    (options, args) = parser.parse_args()
    
    inputArg = options.input
    workingDir = "/global/project/projectdirs/hacc/PDACS/working/"
    if not os.path.exists(workingDir):
        os.makedirs(workingDir)
        os.umask(002)
    #get the name of the file - used to make the temp output files.
    workingName = inputArg.split("/")
    workingName = workingName[len(workingName)-1]

    #print options.outFile
 
    outputFile = workingDir + workingName 
    #workingDir + workingName + ".0"
    print outputFile
    if not os.path.exists(workingDir):
        os.makedirs(workingDir)
    
    flags = ' '
    if options.noRank == 'true':
        flags = ' --no-rank-info '

    pbsPath = create_pbs_file(workingDir, flags, inputArg, outputFile, workingName)

    #use the helper to submit the job
    res, con = NewtHelper.authenticate()
    cookieStr = res['set-cookie']
    res, con = NewtHelper.execute_job(pbsPath, cookieStr)    
    print 'Submitted job: ' + NewtHelper.getJobID(con)
    jobid = NewtHelper.getJobID(con)
    print "Job queued successfully."
    status = NewtHelper.getJobStatus(jobid, cookieStr)
    status = NewtHelper.waitToFinish(jobid, cookieStr)
    
    if status == "C":
        print "Job completed."
    
    #cmd = "tr ' ' '\t' < %s > %s" %(outputFile, options.outFile)
    #cmd = "tr ' ' \\t < %s > %s" %(outputFile, options.outFile)
    cmd = "sed 's/^[ \t]*//' %s | tr -s ' ' '\t' > %s" %(outputFile, options.outFile)
    os.system(cmd)

    #shutil.copyfile(outputFile, options.outFile)

if __name__ == '__main__':
    run()
