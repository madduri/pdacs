import sys
import os
from subprocess import call
import optparse
import httplib2
import urllib, logging
import string
import NewtHelper
from time import time, sleep
import json
import shutil
import sqlite3 as lite

def set_param_file_values(paramFile):
    global numParticles
    global boxSize
    global hubbleConst
    global omegaDM
    
    pfile = open(paramFile, "r")
    for paramLine in pfile.readlines():
        if "=" in paramLine:
            temp = paramLine.split("=")
            val = temp[1].strip()
            if temp[0] == 'hubble':
                hubbleConst = val
            elif temp[0] == 'np':
                numParticles = val
            elif temp[0] == 'box_size':
                boxSize = val
            elif temp[0] == 'Omega_m':
                omegaDM = val

def set_params(paramFile):
    global numParticles
    global boxSize
    global hubbleConst
    global omegaDM
    global n_s
    global w_de
    global sigma_8
    global snapshotname
    global red_shift
    global omega_b  

    con = lite.connect(paramFile)
    cur = con.cursor()
    cur.execute("select value from metadata where name='Snapshot'")
    con.commit()
    row = cur.fetchone()
    snapshotname = row[0]

    con = lite.connect(paramFile)
    cur = con.cursor()
    cur.execute("select value from metadata where name='a_val'")
    con.commit()
    row = cur.fetchone()
    red_shift = row[0]

    con = lite.connect(paramFile)
    cur = con.cursor()
    cur.execute("select value from metadata where name='hubble'")
    con.commit()
    row = cur.fetchone()
    hubbleConst = row[0]
    
    cur.execute("select value from metadata where name='np'")
    con.commit()
    row = cur.fetchone()
    numParticles = row[0]

    cur.execute("select value from metadata where name='box_size [Mpc/h]'")
    con.commit()
    row = cur.fetchone()
    boxSize = row[0]

    cur.execute("select value from metadata where name='Omega_m'")
    con.commit()
    row = cur.fetchone()
    omegaDM = row[0]

    cur.execute("select value from metadata where name='n_s'")
    con.commit()
    row = cur.fetchone()
    n_s = row[0]

    cur.execute("select value from metadata where name='w_de'")
    con.commit()
    row = cur.fetchone()
    w_de = row[0]

    cur.execute("select value from metadata where name='Sigma_8'")
    con.commit()
    row = cur.fetchone()
    sigma_8 = row[0]

    cur.execute("select value from metadata where name='Omega_bar'")
    con.commit()
    row = cur.fetchone()
    omega_b = row[0]

def create_in_file(filePath, inputArg, inputType, distributeType, massFactor, 
                   overloadZoneSize, linkingLength, minParticlesPerHalo, outputParticleFraction, 
                   minNumberParticles, useMCPCenterFinder, useMBPCenterFinder, 
                   outputFOFPorperties, outputSODProperties, outputParticlesPerHalo, 
                   outputParticles, outputHaloCatalog, outputAllParticles, outputAllParticleTags):
    
    inFile = open(filePath, 'w')
    
    inputText = """################################################################################
# Header version information
################################################################################
HALOFINDER_HEADER_VERSION 1.0.0

################################################################################
# Input base name ending in '.' if followed by processor id
################################################################################
INPUT_BASE_NAME %s
""" % inputArg

    inputText = inputText + """
################################################################################
# Input data style (RECORD = .cosmo)  (BLOCK = .gadget2)
################################################################################
INPUT_TYPE %s
""" % inputType

    inputText = inputText + """
################################################################################
# Particle distribution style
#  ROUND_ROBIN indicates particles must be looked at by all processors
#  ONE_TO_ONE  indicates that particles physically reside on matchin processor
################################################################################
DISTRIBUTE_TYPE %s
""" % distributeType

    inputText = inputText + """
################################################################################
# Output base name
################################################################################
OUTPUT_BASE_NAME %s
""" % filePath

    inputText = inputText + """
################################################################################
# Mass factor (all masses read are multiplied by this)
# Default unit is Msun/h
# Any units may be used, but RHO_C must match those units
################################################################################
MASS_CONVERT_FACTOR %s
""" % massFactor

    inputText = inputText + """
################################################################################
# Distance factor (all positions read are multiplied by this)
# Default unit is Mpc/h
# Any units may be used, but RHO_C must match those units
################################################################################
DIST_CONVERT_FACTOR 1.0
    
################################################################################
# RHO_C convert factor
# Default value of 2.77e11 in units of Msun/h and Mpc/h
# Any units may be used, but RHO_C must match those units
################################################################################
RHOC_CONVERT_FACTOR 1.0
    
################################################################################
# SOD_MASS convert factor
# SOD_MASS is used to divide FOF halo mass to get initial SOD radius
# Default value of 1.0e14 in units of Msun/h
# Any units may be used, but SOD_MASS must match those units
################################################################################
SOD_MASS_CONVERT_FACTOR 1.0
    
################################################################################
# Box size (rL)
################################################################################
BOX_SIZE %s
""" % boxSize

    inputText = inputText + """
################################################################################
# Overload zone size (dead zone)
################################################################################
OVERLOAD_SIZE %s
""" % overloadZoneSize

    inputText = inputText + """
################################################################################
# Number of particles in all files (np^3)
################################################################################
NUMBER_OF_PARTICLES %s
""" % numParticles

    inputText = inputText + """
################################################################################
# Minimum distance between particles in a halo (bb)
################################################################################
MINIMUM_PARTICLE_DISTANCE %s
""" % linkingLength

    inputText = inputText + """
################################################################################
# Minimum number of particles in a halo (pmin)
################################################################################
MINIMUM_PARTICLES_PER_HALO %s
""" % minNumberParticles

    inputText = inputText + """
################################################################################
# Omega dm
################################################################################
OMEGADM %s
""" % omegaDM

    inputText = inputText + """
################################################################################
# Hubble constant
################################################################################
HUBBLE_CONSTANT %s
""" % hubbleConst

    inputText = inputText + """
################################################################################
# Subhalo number of particle neighbors used in SPH density calculation
# This number of closest neighbors contributes to the local density
################################################################################
NUM_SPH_DENSITY 64
    
################################################################################
# Subhalo number of particle neighbors used in subhalo grouping
# This number of close neighbors are considered when deciding where
# to place a particle in the candidates list
################################################################################
NUM_SUBHALO_NEIGHBORS 20
    
################################################################################
# Subhalo minimum subhalo size
# Smallest subhalo recognized as a separate unit
################################################################################
MIN_SUBHALO_SIZE 20
    
################################################################################
# Subhalo minimum FOF halo size to have subhalo finder run on
# Subhalo finder won't run on FOF halos smaller than this number
################################################################################
MIN_FOF_SUBHALO 5000
    
################################################################################
# Subhalo alpha factor which controls the cut or grow of subhalos
# at a saddlepoint.  If set to 1.0 always cut the smaller candidate.
# If higher small subhalos and tails grow more aggressively
################################################################################
ALPHA_SUBHALO 1.0
    
################################################################################
# Subhalo beta factor which controls the Poisson noise significance test
# If set to 0.0 then all candidates are deemed significant, if set bigger
# small candidates will be absorbed into large one during subgrouping
################################################################################
BETA_SUBHALO 0.0
    
################################################################################
# Run the MCP (most connected particle) FOF halo center finder
################################################################################
USE_MCP_CENTER_FINDER %s
""" % useMCPCenterFinder

    inputText = inputText + """
################################################################################
# Run the MBP (most bound particle) FOF halo center finder
################################################################################
USE_MBP_CENTER_FINDER %s
""" % useMBPCenterFinder
    
    inputText = inputText + """
################################################################################
# Use the minimum potential array already supplied by the simulation
################################################################################
USE_MINIMUM_POTENTIAL 0
   
################################################################################
# Output all particle data with mass field replaced by halo tag
################################################################################
OUTPUT_PARTICLES %s
""" % outputParticles

    inputText = inputText + """
################################################################################
# Output the halo catalog of one entry per halo (.cosmo and ascii format)
################################################################################
OUTPUT_HALO_CATALOG %s
""" % outputHaloCatalog
    
    inputText = inputText + """
################################################################################
# Output FOF halo properties report (ascii)
################################################################################
OUTPUT_FOF_PROPERTIES %s
""" % outputFOFPorperties
    
    inputText = inputText + """
################################################################################
# Output SOD halo properties report (ascii)
################################################################################
OUTPUT_SOD_PROPERTIES %s
""" % outputSODProperties

    inputText = inputText + """
################################################################################
# Output all particles in halos with at least this number of particles
################################################################################
MINIMUM_PARTICLES_PER_OUTPUT_HALO %s
""" % outputParticlesPerHalo

    inputText = inputText + """
################################################################################
# The fraction of all particles in halos meeting the halo size cut to output
################################################################################
OUTPUT_PARTICLE_FRACTION %s
""" % outputParticleFraction

    inputText = inputText + """
################################################################################
# The minimum number of particles per halo in the fractional output
################################################################################
OUTPUT_MINIMUM_PARTICLES_PER_HALO %s
""" % minParticlesPerHalo

    inputText = inputText + """
################################################################################
# Output all particle tags
################################################################################
OUTPUT_ALL_PARTICLE_TAGS %s
""" % outputAllParticleTags

    inputText = inputText + """
################################################################################
# Output all particles
################################################################################
OUTPUT_ALL_PARTICLES 0
    
################################################################################
# Output Subhalo properties report (ascii)
################################################################################
OUTPUT_SUBHALO_PROPERTIES 0
    
    
"""
        
    inFile.write(inputText)
    inFile.close()
    print 'Created input file.'
    
def create_pbs_file(filePath, nodes, procs, exetime, exeQueue, workingDir, workingName):
    numnodes = int(nodes)
    numprocs = int(procs)
    np = numnodes*numprocs
    exeTime = int(exetime)
    pbsText = """
#PBS -q %(exeQueue)s
#PBS -A hacc
#PBS -l nodes=%(numnodes)d:ppn=%(numprocs)d
#PBS -l walltime=00:%(exeTime)d:00
#PBS -N halotest
#PBS -e %(workingDir)s/halotest.$PBS_JOBID.err
#PBS -o %(workingDir)s/halotest.$PBS_JOBID.out
#PBS -V
    
module unload pgi
module load gcc


cd %(workingDir)s

mpirun -n %(np)d /global/project/projectdirs/hacc/PDACS/HaloFind_v3/HaloTestP %(inputfile)s
chmod 777 %(inputfile)s.*
chgrp hacc %(inputfile)s.*
""" % {"exeQueue":exeQueue, "numnodes":numnodes, "numprocs":numprocs, "exeTime":exeTime, "workingDir":workingDir, "workingDir":workingDir, "workingDir":workingDir, "np":numnodes*numprocs, "inputfile":filePath, "inputfile":filePath, "inputfile":filePath}    
    pbsPath = workingDir + "haloTest." + workingName +".pbs"
    pbsFile = open(pbsPath, 'w')
    pbsFile.write(pbsText)
    pbsFile.close()
    print 'Created PBS file.'
    return pbsPath
    

def run():    
    parser = optparse.OptionParser()
    parser.add_option('--input', action="store", dest="input")
    parser.add_option('--input-type', action="store", dest="inputType")
    parser.add_option('--distribute-type', action="store", dest="distributeType")
    parser.add_option('--massFactor', action="store", dest="massFactor")
    parser.add_option('--overloadZoneSize', action="store", dest="overloadZoneSize")
    parser.add_option('--linkingLength', action="store", dest="linkingLength")
    parser.add_option('--minParticlesPerHalo', action="store", dest="minParticlesPerHalo")
    parser.add_option('--centerFinder', action="store", dest="centerFinder")
    parser.add_option('--minNumberParticles', action="store", dest="minNumberParticles")
    parser.add_option('--numNodes', action="store", dest="numNodes")
    parser.add_option('--numProcs', action="store", dest="numProcs")    
    parser.add_option('--exeTime', action="store", dest="exeTime")    
    parser.add_option('--exeQueue', action="store", dest="exeQueue")    
    parser.add_option('--outputParticleFraction', action="store", dest="outputParticleFraction")
#    parser.add_option('--outParticlesFile', action="store", dest="outParticlesFile")
#    parser.add_option('--outputHaloCatalogFile', action="store", dest="outputHaloCatalogFile")
    parser.add_option('--outputFOFPorpertiesFile', action="store", dest="outputFOFPorpertiesFile")
    parser.add_option('--outputSODPropertiesFile', action="store", dest="outputSODPropertiesFile")
    parser.add_option('--outputParticlesPerHaloFile', action="store", dest="outputParticlesPerHaloFile")
    parser.add_option('--outputAllParticleTagsFile', action="store", dest="outputAllParticleTagsFile")
    parser.add_option('--outputSODPropertyBinsFile', action="store", dest="outputSODPropertyBinsFile")
    parser.add_option('--dbiFile', action="store", dest="dbiFile")
    (options, args) = parser.parse_args()
    
    inputType = options.inputType
    nodes = options.numNodes
    procs = options.numProcs
    distributeType = options.distributeType
    massFactor = options.massFactor
    overloadZoneSize = options.overloadZoneSize
    linkingLength = options.linkingLength
    minParticlesPerHalo = options.minParticlesPerHalo
    outputParticleFraction = options.outputParticleFraction
    minNumberParticles = options.minNumberParticles
#    ParticlesFile = options.outParticlesFile
#    HaloCatalogFile = options.outputHaloCatalogFile
    FOFPropertiesFile = options.outputFOFPorpertiesFile
    SODPropertiesFile = options.outputSODPropertiesFile
    ParticlesPerHaloFile = options.outputParticlesPerHaloFile
    AllParticleTagsFile = options.outputAllParticleTagsFile
    SODPropertyBinsFile = options.outputSODPropertyBinsFile

    useMCPCenterFinder = "0"
    useMBPCenterFinder = "0"

    if options.centerFinder == "MCP":
        useMCPCenterFinder = "1"
    if options.centerFinder == "MBP":
        useMBPCenterFinder = "1"    

    #force all output files to be generated
    outputParticles = "1"
    outputHaloCatalog = "0"
    outputFOFPorperties = "1"
    outputSODProperties = "1"
    outputParticlesPerHalo = "1"
    outputAllParticleTags = "1"
    outputAllParticles = "1"

    workingDir = "/global/project/projectdirs/hacc/PDACS/working/"
    workingName = None

    paramFile = options.input    
    #parse the input file from the directory the snapshot was taken from
    set_params(paramFile)

    inputArg = "/global/project/projectdirs/hacc/PDACS/Coyote/Grid/" + snapshotname + "."
    paramFile = options.input
    workingName = options.dbiFile.split("/")
    workingName = workingName[len(workingName)-1]
    #set the path to create files
    if not os.path.exists(workingDir):
        os.makedirs(workingDir)
    filePath = workingDir + "haloTest."+ workingName + ".in"
    print filePath
    
    #write the values out to an input file
    create_in_file(filePath, inputArg, inputType, distributeType, massFactor, 
                   overloadZoneSize, linkingLength, minParticlesPerHalo, outputParticleFraction, 
                   minNumberParticles, useMCPCenterFinder, useMBPCenterFinder, 
                   outputFOFPorperties, outputSODProperties, outputParticlesPerHalo, 
                   outputParticles, outputHaloCatalog, outputAllParticles, outputAllParticleTags)
    
    #write the pbs file to execute on Carver
    pbsPath = create_pbs_file(filePath, nodes, procs, options.exeTime, options.exeQueue, workingDir, workingName)
    
    #authenticate with NEWT
    res, con = NewtHelper.authenticate()
    cookieStr = res['set-cookie']
    print 'Authenticated'  
    print cookieStr
    res, con = NewtHelper.execute_job(pbsPath, cookieStr)
    print 'Submitted job: ' + NewtHelper.getJobID(con)
    jobid = NewtHelper.getJobID(con)
    print "Job queued successfully."
    status = NewtHelper.getJobStatus(jobid, cookieStr)
    status = NewtHelper.waitToFinish(jobid, cookieStr)
    if status == "C":
        print "Job completed."
    
    #output the files
    #TODO what is this file extension when made for these two?
    #shutil.copyfile(filePath + ".bighaloparticles",options.outParticlesFile)
    shutil.copyfile(filePath + ".haloparticles",options.outputParticlesPerHaloFile)
    shutil.copyfile(filePath + ".fofproperties",options.outputFOFPorpertiesFile)
    shutil.copyfile(filePath + ".sodproperties",options.outputSODPropertiesFile)
    #shutil.copyfile(filePath + ".halocatalog",options.outputHaloCatalogFile)
    shutil.copyfile(filePath + ".haloparticletags",options.outputAllParticleTagsFile)
    shutil.copyfile(filePath + ".sodpropertybins",options.outputSODPropertyBinsFile)

    os.remove(filePath + ".haloparticles")
    os.remove(filePath + ".fofproperties")
    os.remove(filePath + ".sodproperties")
    os.remove(filePath + ".haloparticletags")
    os.remove(filePath + ".sodpropertybins")
    os.remove(filePath + ".bighaloparticles")
    os.remove(filePath)
    os.remove(pbsPath)

    dbifilename = options.dbiFile
    
    fullpath = snapshotname.split("/")
    model = fullpath[0]
    size = fullpath[1]
    realization = fullpath[2]
    snapshot = fullpath[3]

    redshift = 1/red_shift - 1
#    print red_shift
    cambFile = "/global/project/projectdirs/hacc/PDACS/Coyote/Grid/" + model + "/" + size + "/" + realization + "/" + "camb" + model + ".tf" 
    con = lite.connect(options.dbiFile)
    cur = con.cursor() 

    cur.execute("create table metadata (name TEXT, value NUMERIC)")
    cur.execute("insert into metadata values( ?, ?)", ( 'program', 'HaloFinder' ))
    cur.execute("insert into metadata values( ?, ? )", ( 'Simulation', 'Coyote' ))
    cur.execute("insert into metadata values( ?, ? )", ( 'Type', 'Grid' ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'Model', model ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'Size [Mpc]', size ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'Realization', realization))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'Snapshot', snapshot ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'hubble', hubbleConst ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'box_size [Mpc/h]', boxSize ))
#    cur.execute("insert into metadata values( ?, ? ) ", ( 'seed', boxSize ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'numParticles', numParticles))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'Omega_m', omegaDM))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'Omega_b', omega_b))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'ns', n_s))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'w0', w_de))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'sigma8', sigma_8))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'InputType', inputType ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'Nodes', nodes ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'Procs', procs ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'distributeType', distributeType ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'massFactor', massFactor ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'overloadZoneSize', overloadZoneSize ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'linkingLength', linkingLength ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'minParticlesPerHalo', minParticlesPerHalo ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'red_shift', redshift ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'outputParticleFraction', outputParticleFraction ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'minNumberParticles', minNumberParticles ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'camboutput', cambFile))
    #cur.execute("insert into metadata values( ?, ? ) ", ( 'ParticlesFile', ParticlesFile))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'FOFPropertiesFile', FOFPropertiesFile ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'SODPropertiesFile', SODPropertiesFile ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'ParticlesPerHaloFile', ParticlesPerHaloFile ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'AllParticleTagsFile', AllParticleTagsFile ))
    cur.execute("insert into metadata values( ?, ? ) ", ( 'SODPropertyBinsFile', SODPropertyBinsFile ))
    con.commit()

if __name__ == '__main__':
    run()
