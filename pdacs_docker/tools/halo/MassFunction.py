#!/usr/bin/python

import sqlite3 as lite
import sys, os, os.path

#some useful piece of code to make sure that MF, cM, fsigma, growth_ode are considered as python modules from whereever they exist in the code
#PYTHON_TOOLS_PATH is an environment variable needed to let the "wrapper" know where rest of the modules or the scripts are. 
#If there are multiple files/modules required by the main python script, they should have __init__.py in the same folder so that wrapper or main tool
#could look it up. 

#sys.path.append(os.path.join(os.path.dirname(__file__), os.environ["PYTHON_TOOLS_PATH"]))
sys.path.append(os.path.join(os.path.dirname(__file__), '/global/project/projectdirs/hacc/PDACS/Halo_analysis_v2'))

from MF import *
from cM import *
from fsigma import *
from growth_ode import *
from optparse import OptionParser
import time

#parsing input arguments
parser = OptionParser()
parser.add_option("-j", "--fType", dest="fType", type="int", default=1, action="store")
parser.add_option("-m", "--mBins", dest="massbins", type="int", default=20, action="store", help="Number of Mass bins[%default]")
parser.add_option("-n", "--nFiles", dest="numFiles", type="int", default=1, action="store", help="Number of Mass bins[%default]")
parser.add_option("-z", "--zSeek", dest="zseek", type="float", default=0.0, action="store", help="Number of Mass bins[%default]")
parser.add_option("-l", "--mNum", dest="model_num", type="int", default=0, action="store", help="Number of Mass bins[%default]")
parser.add_option("-e", "--hubble", dest="hubble", type="float", default=0.5977, action="store", help="Number of Mass bins[%default]")
parser.add_option("-g", "--omegaM", dest="Omegam", type="float", default=0.431, action="store", help="Number of Mass bins[%default]")
parser.add_option("-a", "--ns", dest="ns", type="float", default=0.9468, action="store", help="Number of Mass bins[%default]")
parser.add_option("-s", "--sigma8", dest="sigma8", type="float", default=0.8161, action="store", help="Number of Mass bins[%default]")
parser.add_option("-w", "--w0", dest="w0", type="float", default=-0.816, action="store", help="Number of Mass bins[%default]")
parser.add_option("-b", "--boxsize", dest="boxsize", type="float", default=107.6, action="store", help="Number of Mass bins[%default]")
parser.add_option("-p", "--particlenum", dest="particlenum", type="int", default=512, action="store", help="Number of Mass bins[%default]")
parser.add_option("-r", "--minparticle", dest="minparticle", type="int", default=1000, action="store", help="Number of Mass bins[%default]")
parser.add_option("-d", "--delta", dest="Delta", type="float", default=200.0, action="store", help="Number of Mass bins[%default]")
parser.add_option("-u", "--haloFile", dest="halofile", type="string", default="halo_m001.gadget", action="store", help="Tk File[%default]")
parser.add_option("-i", "--inputFile", action="store", dest="inputFile", type="string")
parser.add_option("-t", "--tkFile", dest="tkFilename", type="string", default="cambM001.tf", action="store", help="Tk File[%default]")
parser.add_option("-c", "--dirPath", dest="dirPathname", type="string", default="/global/project/projectdirs/hacc/PDACS/Halo_analysis_v2", action="store", help="Input dir Path [%default]")
parser.add_option("-o", "--output", dest="outFilename", action="store", help="SQLite db name [%default]")

(options, args) = parser.parse_args()
#parser.print_help()

#if len(args) != 2:
#	parser.print_help()
#end of parsing arguments

rhoc = 2.77536627e11; #Msun.h^2/Mpc^3

tin = time.time()

#print os.getcwd()
dirPname = options.dirPathname
print options.outFilename
#infile = dirPname+'/'+options.inFilename  
#inputparams = read_input_params(infile) 

Tkfile= dirPname + '/' + options.tkFilename 

massbins= options.massbins
numFiles= options.numFiles
zseek= options.zseek
model_num= options.model_num
hubble= options.hubble
Omegam= options.Omegam
ns= options.ns
sigma8= options.sigma8
w0= options.w0
print "hubble = ", hubble
print "OmegaM = ", Omegam
print "ns = ", ns
print "sigma8 = ", sigma8
print "w0 = ", w0
boxsize= options.boxsize
particlenum= options.particlenum
minparticle= options.minparticle
Delta= options.Delta
#these three inputs should be fed from command line as well?
fileroot= dirPname
halofile= options.halofile

massres= Omegam*rhoc*(boxsize/particlenum)**3 #Msun/h
simvol= (boxsize)**3 #Mpc^3/h^3
print 'mass res= %le' % (massres), 'Msun/h'    
    
k,Tk= read_Tk(Tkfile)
N= Pk_norm(sigma8, ns, k, Tk, hubble)
scalefac, D0, D1= Da_ode(Omegam, 1-Omegam, w0)
Ds, lDs= interp_D(zseek,scalefac,D0, D1)
print "D(z) and log D/log z(z)= ", Ds, lDs, 'at z= ',zseek
print "time takes= ",time.time()-tin

#connect to the db 
con = lite.connect(options.outFilename)
cur = con.cursor()

if options.fType == 1:
        fofproperties= "fofproperties"
        halofileroot=fileroot+'/'+halofile
        fofcount_col= readheaders(halofileroot, fofproperties,'fof_halo_count')
        FOFcount= readhalofiles(halofileroot, fofproperties, fofcount_col, numFiles)
         
        print "\ntotal # of fof clusters= ",len(FOFcount)
        FOFbinstart=min(FOFcount)*massres
        FOFbinend=max(FOFcount)*massres
        FOFmass=FOFcount*(1.0-1.0/FOFcount**0.65)*massres
        print "minFOFcount  maxFOFcount" ,min(FOFcount), max(FOFcount)
	print "mass Msun/h= ", '%le %le' %(FOFbinstart, FOFbinend)
	cur.execute("drop table if exists fof_mf") 
	cur.execute("CREATE TABLE fof_mf(FOF_Mass_Msun_by_h REAL,\n num_clusters INTEGER,\n dn_by_dln REAL,\n frac_err REAL,\n One_by_sigma_M REAL,\n f_sigma REAL,\n fsigma_fit REAL)")        
	meanmass, hist, binsize, bin_edges=calc_mf(FOFbinstart, FOFbinend, massbins, FOFmass)
	for i in range(0,massbins-1):  
           sigmaM= sigmam(Omegam, ns, meanmass[i], N, k, Tk, hubble)
           logsigmaM= logsigm(Omegam, ns, meanmass[i], N, k, Tk, sigmaM, hubble)
           sigmaM= Ds*sigmaM
	   tuples= (meanmass[i], int(hist[i]), hist[i]/binsize/simvol,1/hist[i]**0.5, 1/sigmaM, hist[i]/binsize/simvol*meanmass[i]/(-rhoc*Omegam*logsigmaM), MF_fit(sigmaM, zseek))
	   #print tuples 
           if hist[i] >0: cur.execute("insert into fof_mf values(?, ?, ?, ?, ?, ?, ?)", tuples)
	con.commit()  

if options.fType == 2 or options.fType == 3:
        tin= time.time()
        sodproperties= "sodproperties"
        halofileroot=fileroot+'/'+halofile
        socount_col= readheaders(halofileroot, sodproperties,'sod_halo_count')
        SOcount= readhalofiles(halofileroot, sodproperties, socount_col, numFiles)
        print "\ntotal # of SO clusters= ",len(SOcount)
        SObinstart=max(min(SOcount), minparticle)*massres
        SObinend=max(SOcount)*massres
        SOmass=SOcount*massres
        SOradius= (3.0*SOmass/4.0/pi/rhoc/Delta)**(1.0/3.0)
        print "minSO maxSO"
        print "# particles= ",min(SOcount), max(SOcount)
        print "mass Msun/h= ", '%le %le' %(SObinstart, SObinend)
        cur.execute("drop table if exists sod_mf")
        cur.execute("CREATE TABLE sod_mf(SO_Mass_Msun_h REAL,\n num_clusters INTEGER,\n dn_dln_M REAL,\n frac_err REAL,\n One_by_sigma_M REAL,\n f_sigma REAL)")
        meanmass, hist, binsize, bin_edges= calc_mf(SObinstart, SObinend, massbins, SOmass)
        for i in range(0,massbins-1):
           sigmaM= sigmam(Omegam, ns, meanmass[i], N, k, Tk, hubble)
           logsigmaM= logsigm(Omegam, ns, meanmass[i], N, k, Tk, sigmaM, hubble)
           sigmaM= Ds*sigmaM
           if hist[i] >2:
                sod_tuples = (meanmass[i], int(hist[i]), hist[i]/binsize/simvol,1/hist[i]**0.5, 1/sigmaM, hist[i]/binsize/simvol*meanmass[i]/(-rhoc*Omegam*logsigmaM))
                cur.execute("insert into sod_mf values(?, ?, ?, ?, ?, ?)", sod_tuples)
        print "time takes= ",time.time()-tin
        con.commit()

if options.fType == 3:
        def fitz(z, var):
          if var ==0:
            if z<=1: y= 1.08*(1+z)**0.27
            if z>1: y= 1.3*z**-0.19
          if var>0 and z<=1 and z>0: y=1.08*(1+z)**0.35
          if var>0 and z==0: y=0.95
          return y
        
	tin= time.time()
        sodprofile= "sodpropertybins"
        halofileroot=fileroot+'/'+halofile
        cts_col= readheaders(halofileroot, sodprofile,'sod_halo_bin_count')
        radius_col= readheaders(halofileroot, sodprofile,'sod_halo_bin_radius')
        overden_col= readheaders(halofileroot, sodprofile,'sod_halo_bin_rho_ratio')
        print "\ndoing c-M calculation now.."
        conc, concerr, norm, massSO= conc_each_halo_lessmem(SOmass, SOradius, halofileroot, sodprofile, numFiles, SObinstart, cts_col, radius_col, overden_col)
        cmean, cmeanerr, meanmass, count, variance= get_bin_cmean(massSO, conc, concerr, bin_edges, massbins)
        cur.execute("drop table if exists cM")
        cur.execute("CREATE TABLE cM(SO_Mass_Msun_h REAL,\n num_clusters INTEGER,\n cmean REAL,\n frac_err REAL,\n variance REAL)")
        for i in range(0,massbins-1):
            if (count[i] >3 and variance[i] !=0.0 and cmean[i] !=0):
                cm_tuple = (meanmass[i], int(count[i]),cmean[i]*fitz(zseek, model_num), (cmeanerr[i]**2+cmean[i]**2/count[i])**0.5, (variance[i]/cmean[i]**2-1.0)**0.5)
                cur.execute("insert into cM values(?, ?, ?, ?, ?)", cm_tuple)
        print "time takes= ",time.time()-tin
        con.commit()

