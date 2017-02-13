import sys
import os
import optparse
import shutil
import NewtHelper

def create_pbs_file(toolDir, workingDir, m, b, ns, w, s, z, dbFile):
    pbsText = """
#PBS -q debug
#PBS -A hacc
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:15:00
#PBS -N cmTool
#PBS -e cmTool.$PBS_JOBID.err
#PBS -o cmTool.$PBS_JOBID.out
#PBS -V

module unload pgi
module load gcc
module load gsl

cd %s

/global/project/projectdirs/hacc/PDACS/cm/cM.sh %s %s %s %s %s %s %s

""" % (toolDir, m, b, ns, w, s, z, dbFile)
    pbsPath = workingDir + "cmTool.pbs"
    pbsFile = open(pbsPath, 'w')
    pbsFile.write(pbsText)
    pbsFile.close()
    print 'Created PBS file.'
    return pbsPath


def run():
    parser = optparse.OptionParser()
    parser.add_option('--m', action="store", dest="m")
    parser.add_option('--b', action="store", dest="b")
    parser.add_option('--ns', action="store", dest="ns")
    parser.add_option('--w', action="store", dest="w")
    parser.add_option('--s', action="store", dest="s")
    parser.add_option('--z', action="store", dest="z")
    parser.add_option('--outFile', action="store", dest="outFile")
    (options, args) = parser.parse_args()

    toolDir = "/global/project/projectdirs/hacc/PDACS/cm/"
    if not os.path.exists(toolDir):
        os.makedirs(toolDir)

    workingDir = "/global/project/projectdirs/hacc/PDACS/working/"
    if not os.path.exists(workingDir):
        os.makedirs(workingDir)

    #get the name of the file - used to make the temp output files.
    workingName = options.outFile.split("/")
    workingName = workingName[len(workingName)-1]

    outputFile = workingDir + workingName + ".sqlite"
    if not os.path.exists(workingDir):
        os.makedirs(workingDir)

    pbsPath = create_pbs_file(toolDir, workingDir, options.m, options.b, options.ns, options.w, options.s, options.z, outputFile)
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

    os.chmod(outputFile, 0777)
    #os.chmod(options.outFile, 0777)
    print outputFile
    print options.outFile
    shutil.copyfile(outputFile, options.outFile)
    print 'Copied file to output.'

if __name__ == '__main__':
    run()

