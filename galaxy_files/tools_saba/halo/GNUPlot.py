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

def create_plt_file(workingDir, workingName, inputArg):    
    pltText = """
    set term postscript
    set output '%s.ps'
    plot"%s" u 3:($2*1e10) w p
    
    """ % (workingName, inputArg)
    
    pltPath = workingDir + workingName + ".plt"
    pltFile = open(pltPath, 'w')
    pltFile.write(pltText)
    pltFile.close()
    return pltPath

#def create_pbs_file(workingDir, pltPath, workingName):
#    pbsText = """
#    #PBS -q debug
#    #PBS -A hacc
#    #PBS -l nodes=1:ppn=1
#    #PBS -l walltime=00:15:00
#    #PBS -N gnuplot
#    #PBS -e gnuplot.$PBS_JOBID.err
#    #PBS -o gnuplot.$PBS_JOBID.out
#    #PBS -V
#    
#    cd %s
#    
#    module load gnuplot
#    
##    gnuplot %s
#    
#    """ % (workingDir, pltPath)
#    pbsPath = workingDir + workingName + ".pbs"
#    pbsFile = open(pbsPath, 'w')
#    pbsFile.write(pbsText)
#    pbsFile.close()
#    print 'Created PBS file.'
#    return pbsPath


def run():

    parser = optparse.OptionParser()
    parser.add_option('--input', action="store", dest="input")
    parser.add_option('--outFile', action="store", dest="outFile")
    (options, args) = parser.parse_args()
    
    inputArg = options.input
    
    #workingDir = "/project/projectdirs/hacc/PDACS/working/"
    workingDir = "./"    

    #if not os.path.exists(workingDir):
    #    os.makedirs(workingDir)
    
    #get the name of the file - used to make the temp output files.
    workingName = inputArg.split("/")
    workingName = workingName[len(workingName)-1]
    
    filePath = workingDir + workingName + ".ascii"
    outputFile = workingDir + workingName + ".ps"

    #create the required files
    pltPath = create_plt_file(workingDir, workingName, inputArg)
    
    #pbsPath = create_pbs_file(workingDir, pltPath, workingName)
    
    #use the helper to submit the job
    #res, con = NewtHelper.authenticate()
    #cookieStr = res['set-cookie']
    #res, con = NewtHelper.execute_job(pbsPath, cookieStr)    
    #print 'Submitted job: ' + NewtHelper.getJobID(con)
    #jobid = NewtHelper.getJobID(con)
    #print "Job queued successfully."
    #status = NewtHelper.getJobStatus(jobid, cookieStr)
    #status = NewtHelper.waitToFinish(jobid, cookieStr)
    
    #if status == "C":
    #    print "Job completed."
    #
    #shutil.copyfile(outputFile, options.outFile)
    cmd = "module load gnuplot; gnuplot %s" % pltPath

if __name__ == '__main__':
    run()

