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

def create_in_file(filePath, inputArg, inputType, distributeType, massFactor, 
                   overloadZoneSize, linkingLength, minParticlesPerHalo, deut, outputParticleFraction, 
                   minNumberParticles, useMCPCenterFinder, useMBPCenterFinder, outputParticles, 
                   outputHaloCatalog, outputFOFPorperties, outputSODProperties, outputParticlesPerHalo, 
                   outputAllParticles, outputAllParticleTags):
    
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
""" % minParticlesPerHalo

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
# Deut
################################################################################
DEUT %s
""" % deut

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
""" % minNumberParticles

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
    
def create_pbs_file(filePath, nodes, procs, workingDir, workingName):
    pbsText = """
#PBS -q debug
#PBS -A hacc
#PBS -l nodes=%s:ppn=%s
#PBS -l walltime=00:15:00
#PBS -N halotest
#PBS -e halotest.$PBS_JOBID.err
#PBS -o halotest.$PBS_JOBID.out
#PBS -V
    
module unload pgi
module load gcc
    
cd %s
    
mpirun -np 32 /global/project/projectdirs/hacc/PDACS/HaloFind_v3/HaloTestP %s
""" % (nodes, procs, workingDir, filePath)
    pbsPath = workingDir + workingName +".pbs"
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
    parser.add_option('--deut', action="store", dest="deut")
    parser.add_option('--centerFinder', action="store", dest="centerFinder")
    parser.add_option('--minNumberParticles', action="store", dest="minNumberParticles")
    parser.add_option('--numNodes', action="store", dest="numNodes")
    parser.add_option('--numProcs', action="store", dest="numProcs")    
    parser.add_option('--outputParticleFraction', action="store", dest="outputParticleFraction")
    #parser.add_option('--outputParticles', action="store", dest="outputParticles")
    #parser.add_option('--outputHaloCatalog', action="store", dest="outputHaloCatalog")
    #parser.add_option('--outputFOFPorperties', action="store", dest="outputFOFPorperties")
    #parser.add_option('--outputSODProperties', action="store", dest="outputSODProperties")
    #parser.add_option('--outputParticlesPerHalo', action="store", dest="outputParticlesPerHalo")
    #parser.add_option('--minNumberParticles', action="store", dest="minNumberParticles")
    #parser.add_option('--outputAllParticleTags', action="store", dest="outputAllParticleTags")
    #parser.add_option('--outputAllParticles', action="store", dest="outputAllParticles")
    #parser.add_option('--numNodes', action="store", dest="numNodes")
    #parser.add_option('--numProcs', action="store", dest="numProcs")    

    parser.add_option('--outParticlesFile', action="store", dest="outParticlesFile")
    parser.add_option('--outputHaloCatalogFile', action="store", dest="outputHaloCatalogFile")
    parser.add_option('--outputFOFPorpertiesFile', action="store", dest="outputFOFPorpertiesFile")
    parser.add_option('--outputSODPropertiesFile', action="store", dest="outputSODPropertiesFile")
    parser.add_option('--outputParticlesPerHaloFile', action="store", dest="outputParticlesPerHaloFile")
    parser.add_option('--outputAllParticleTagsFile', action="store", dest="outputAllParticleTagsFile")
    parser.add_option('--outputSODPropertyBinsFile', action="store", dest="outputSODPropertyBinsFile")

    #parser.add_option('--outputAllParticlesFile', action="store", dest="outputAllParticlesFile")
    (options, args) = parser.parse_args()
    
    inputType = options.inputType
    nodes = options.numNodes
    procs = options.numProcs
    distributeType = options.distributeType
    massFactor = options.massFactor
    overloadZoneSize = options.overloadZoneSize
    linkingLength = options.linkingLength
    minParticlesPerHalo = options.minParticlesPerHalo
    deut = options.deut
    outputParticleFraction = options.outputParticleFraction
    minNumberParticles = options.minNumberParticles

    useMCPCenterFinder = "0"
    useMBPCenterFinder = "0"
    #outputParticles = "0"
    #outputHaloCatalog = "0"
    #outputFOFPorperties = "0"
    #outputSODProperties = "0"
    #outputParticlesPerHalo = "0"
    #outputAllParticleTags = "0"
    #outputAllParticles = "0"

    if options.centerFinder == "MCP":
        useMCPCenterFinder = "1"
    if options.centerFinder == "MBP":
        useMBPCenterFinder = "1"    
    #if options.outputParticles == "true":
    #    outputParticles = "1"
    #if options.outputHaloCatalog == "true":
    #    outputHaloCatalog = "1"
    #if options.outputFOFPorperties == "true":
    #    outputFOFPorperties = "1"
    #if options.outputSODProperties == "true":
    #    outputSODProperties = "1"
    #if options.outputParticlesPerHalo == "true":
    #    outputParticlesPerHalo = "1"
    #if options.outputAllParticles == "true":
    #    outputAllParticles = "1"
    #if options.outputAllParticleTags == "true":
    #    outputAllParticleTags = "1"
   


    #force all output files to be generated
    outputParticles = "1"
    outputHaloCatalog = "1"
    outputFOFPorperties = "1"
    outputSODProperties = "1"
    outputParticlesPerHalo = "1"
    outputAllParticleTags = "1"
    outputAllParticles = "1"


    workingDir = "/global/project/projectdirs/hacc/PDACS/working/"
    workingName = None

    inputArg = options.input
    split = os.path.splitext(options.input)
    print inputArg
    inputArg = split[0] + "."
    paramFile = os.path.dirname(options.input) + "/params.ini"
    print inputArg

    #get the name of the file - used to make the temp output files.
    workingName = inputArg.split("/")
    print workingName
    workingName = workingName[len(workingName)-1]
    print workingName

    #parse the input file from the directory the snapshot was taken from
    set_param_file_values(paramFile)
    
    #set the path to create files
    if not os.path.exists(workingDir):
        os.makedirs(workingDir)
    filePath = workingDir + workingName + "in"

    print filePath
    
    #write the values out to an input file
    create_in_file(filePath, inputArg, inputType, distributeType, massFactor, 
                   overloadZoneSize, linkingLength, minParticlesPerHalo, deut, outputParticleFraction, 
                   minNumberParticles, useMCPCenterFinder, useMBPCenterFinder, outputParticles, 
                   outputHaloCatalog, outputFOFPorperties, outputSODProperties, outputParticlesPerHalo, 
                   outputAllParticles, outputAllParticleTags)
    
    #write the pbs file to execute on Carver
    pbsPath = create_pbs_file(filePath, nodes, procs, workingDir, workingName)
    
    #authenticate with NEWT
    res, con = NewtHelper.authenticate()
    cookieStr = res['set-cookie']
    print 'Authenticated'  
    print res
    print con  
    res, con = NewtHelper.execute_job(pbsPath, cookieStr)
    print res
    print con
    print 'Submitted job: ' + NewtHelper.getJobID(con)
    jobid = NewtHelper.getJobID(con)
    print "Job queued successfully."
    status = NewtHelper.getJobStatus(jobid, cookieStr)
    status = NewtHelper.waitToFinish(jobid, cookieStr)
    
    if status == "C":
        print "Job completed."
    
    #output the files
    #TODO what is this file extension when made for these two?
    shutil.copyfile(filePath + ".bighaloparticles",options.outParticlesFile)
    shutil.copyfile(filePath + ".haloparticles",options.outputParticlesPerHaloFile)

    shutil.copyfile(filePath + ".fofproperties",options.outputFOFPorpertiesFile)
    shutil.copyfile(filePath + ".sodproperties",options.outputSODPropertiesFile)
    shutil.copyfile(filePath + ".halocatalog",options.outputHaloCatalogFile)
    shutil.copyfile(filePath + ".haloparticletags",options.outputAllParticleTagsFile)
    shutil.copyfile(filePath + ".sodpropertybins",options.outputSODPropertyBinsFile)


    #if options.outParticlesFile != "None":
    #    shutil.copyfile(filePath + ".particles",options.outParticlesFile)
    #if options.outputParticlesPerHaloFile != "None":
    #    shutil.copyfile(filePath + ".particlesperhalo",options.outputParticlesPerHaloFile)
    
    #if options.outputFOFPorpertiesFile != "None":
    #    shutil.copyfile(filePath + ".fofproperties",options.outputFOFPorpertiesFile)
    #if options.outputSODPropertiesFile != "None":
    #    shutil.copyfile(filePath + ".sodproperties",options.outputSODPropertiesFile)
    #if options.outputHaloCatalogFile != "None":
    #    shutil.copyfile(filePath + ".halocatalog",options.outputHaloCatalogFile)
    #if options.outputAllParticleTagsFile != "None":
    #    shutil.copyfile(filePath + ".haloparticletags",options.outputAllParticleTagsFile)




if __name__ == '__main__':
    run()
