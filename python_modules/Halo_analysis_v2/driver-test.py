#!/usr/bin/python

import sqlite3 as lite
import sys, os, os.path

#some useful piece of code to make sure that MF, cM, fsigma, growth_ode are considered as python modules from whereever they exist in the code
#PYTHON_TOOLS_PATH is an environment variable needed to let the "wrapper" know where rest of the modules or the scripts are. 
#If there are multiple files/modules required by the main python script, they should have __init__.py in the same folder so that wrapper or main tool
#could look it up. 

#I am modifying this file for some test purposes and cleaning it up, after that people can use this version to create a newer version of this analysis code. 

#sys.path.append(os.path.join(os.path.dirname(__file__), os.environ["PYTHON_TOOLS_PATH"]))
sys.path.append(os.path.join(os.path.dirname(__file__), '/project/projectdirs/hacc/PDACS/Halo_analysis_v2'))

from MF import *
from cM import *
from fsigma import *
from growth_ode import *
from optparse import OptionParser
import time

parser = OptionParser()
parser.add_option("-j", "--fType", dest="fType", type="int", default=1, action="store")
parser.add_option("-c", "--cMFile", dest="cMFile", default="/project/projectdirs/hacc/PDACS/Halo_analysis_v2/halo_m001.gadget.sodpropertybins.0", action="store")
parser.add_option("-m", "--massBins", action="store", type="int", dest="massBins")
parser.add_option("-i", "--inputFile", action="store", dest="inputFile")
parser.add_option("-p", "--propFile", action="store", dest="propFile")
parser.add_option("-o","--outfile", action="store", dest="outfile")
#parser.add_option('--coutfile', action="store", dest="coutfile")
(options, args) = parser.parse_args()

rhoc= 2.77536627e11; #Msun.h^2/Mpc^3

tin= time.time()
#infile= options.inputFile
#print options.inputFile
#inputparams= read_input_params(infile) 
con = lite.connect(options.inputFile)
cur = con.cursor()
#case= int(inputparams[0])
massbins=int(options.massBins)
filenum=1 #there is always one file in current setup and files are so small that there is no reason to split them up that is why it is hardcoded for now to onw
cur.execute("select value from metadata where name='red_shift'")
row = cur.fetchone()
zseek = float(row[0])
model_no=0
cur.execute("select value from metadata where name='hubble'")
row = cur.fetchone()
hubble= float(row[0])
cur.execute("select value from metadata where name='Omega_m'")
row = cur.fetchone()
Omegam= float(row[0])
cur.execute("select value from metadata where name='ns'")
row = cur.fetchone()
ns= float(row[0])
cur.execute("select value from metadata where name='sigma8'")
row = cur.fetchone()
sigma8= float(row[0])
cur.execute("select value from metadata where name='w0'")
row = cur.fetchone()
w0= float(row[0])

#if options.fType == 1:
#   cur.execute("select value from metadata where name='FOFPropertiesFile'")
#   propinFile=cur.fetchone()
#   os.system("/project/projectdirs/hacc/PDACS/pdacs-test/tools/halo/GenericIOPrint.py propinFile propFile 0")

#if options.fType == 2 or options.fType == 3:
#   cur.execute("select value from metadata where name='SODPropertiesFile'")
#   propinFile=cur.fetchone()
#   os.system("/project/projectdirs/hacc/PDACS/pdacs-test/tools/halo/GenericIOPrint.py propinFile propFile 0")
#   cur.execute("select value from metadata where name='SODPropertyBinsFile'")
#   propinFile=cur.fetchone()
#   os.system("/project/projectdirs/hacc/PDACS/pdacs-test/tools/halo/GenericIOPrint.py propinFile cMFile 0")

print hubble, Omegam, ns, sigma8, w0

cur.execute("select value from metadata where name='box_size [Mpc/h]'")
row = cur.fetchone()
boxsize= float(row[0])
cur.execute("select value from metadata where name='numParticles'")
row = cur.fetchone()
particlenum= int(row[0])
cur.execute("select value from metadata where name='minNumberParticles'")
row = cur.fetchone()
minparticle= int(row[0])
print minparticle
Delta= 200
#hard coded for now
#float(inputparams[14])
#fileroot= inputparams[15].rstrip()
#halofile= inputparams[16].rstrip()
cur.execute("select value from metadata where name='camboutput'")
row = cur.fetchone()
#Tkfile="/global/project/projectdirs/hacc/PDACS/Coyote/Grid/M001/L180/G001/cambM001.tf" 
Tkfile=row[0].strip()

massres= Omegam*rhoc*(boxsize/particlenum)**3 #Msun/h
simvol= (boxsize)**3 #Mpc^3/h^3
print 'mass res= %le' % (massres), 'Msun/h'    
print Tkfile
    
k,Tk= read_Tk(Tkfile)
N= Pk_norm(sigma8, ns, k, Tk, hubble)
scalefac, D0, D1= Da_ode(Omegam, 1-Omegam, w0)
Ds, lDs= interp_D(zseek,scalefac,D0, D1)
print "D(z) and log D/log z(z)= ", Ds, lDs, 'at z= ',zseek
print "time takes= ",time.time()-tin
print options.outfile
halofileroot=options.propFile
print halofileroot

if options.fType == 1:
        print "\n################## FOF MF ##########################################\n"
        fofproperties= "fofproperties".rstrip()
        #halofileroot=fileroot+'/'+halofile
        fofcount_col= readheaders(halofileroot, fofproperties,'fof_halo_count')
        FOFcount= readhalofiles(halofileroot, fofproperties, fofcount_col,filenum)
        
        print "\ntotal # of fof clusters= ",len(FOFcount)
        FOFbinstart=min(FOFcount)*massres
        FOFbinend=max(FOFcount)*massres
        FOFmass=FOFcount*(1.0-1.0/FOFcount**0.65)*massres
        print "minFOF maxFOF"
        print "# particles= ",min(FOFcount), max(FOFcount)
	print "mass Msun/h= ", '%le %le' %(FOFbinstart, FOFbinend)
        file= open(options.outfile,'w')
        print >>file, "# FOF Mass Msun/h\t# clusters\tdn/dln M (h/Mpc)^3\tfrac err\t1/sigma(M)\tf(sigma)\tfsigma_fit\t\n"
        meanmass, hist, binsize, bin_edges= calc_mf(FOFbinstart, FOFbinend, massbins, FOFmass)
        for i in range(0,massbins-1):  
           sigmaM= sigmam(Omegam, ns, meanmass[i], N, k, Tk, hubble)
           logsigmaM= logsigm(Omegam, ns, meanmass[i], N, k, Tk, sigmaM, hubble)
           sigmaM= Ds*sigmaM
           if hist[i] >0:  print >> file, '%le\t%5d\t%7.4le\t%7.4lf\t%7.4lf\t%7.4lf\t%7.4lf' %(meanmass[i], hist[i], hist[i]/binsize/simvol,1/hist[i]**0.5, 1/sigmaM, hist[i]/binsize/simvol*meanmass[i]/(-rhoc*Omegam*logsigmaM), MF_fit(sigmaM, zseek))
	file.close()

if options.fType == 2 or options.fType == 3:
        print "\n##################### SO MF ########################################\n"
        tin= time.time()
        sodproperties= "sodproperties".rstrip()
       # halofileroot=fileroot+'/'+halofile
        socount_col= readheaders(halofileroot, sodproperties,'sod_halo_count')
        print socount_col
        SOcount= readhalofiles(halofileroot, sodproperties, socount_col, filenum)   
	print "\ntotal # of SO clusters= ",len(SOcount)
        SObinstart=max(min(SOcount), minparticle)*massres
        SObinend=max(SOcount)*massres
        SOmass=SOcount*massres
        SOradius= (3.0*SOmass/4.0/pi/rhoc/Delta)**(1.0/3.0)
	print "minSO maxSO"
        print "# particles= ",min(SOcount), max(SOcount)
	print "mass Msun/h= ", '%le %le' %(SObinstart, SObinend)
        
        if options.fType == 2: 
	   file= open(options.outfile,'w')
           print>> file, "# SO Mass Msun/h\t# clusters\tdn/dln M(h/Mpc)^3\tfrac err\t1/sigma(M)\tf(sigma)\n"
        meanmass, hist, binsize, bin_edges= calc_mf(SObinstart, SObinend, massbins, SOmass)
        for i in range(0,massbins-1):
           sigmaM= sigmam(Omegam, ns, meanmass[i], N, k, Tk, hubble)
           logsigmaM= logsigm(Omegam, ns, meanmass[i], N, k, Tk, sigmaM, hubble)
	   sigmaM= Ds*sigmaM
           if hist[i] >2 and options.fType == 2:  print>> file, '%le\t%5d\t%7.4le\t%7.4lf\t%7.4lf\t%7.4lf' %(meanmass[i], hist[i], hist[i]/binsize/simvol,1/hist[i]**0.5, 1/sigmaM, hist[i]/binsize/simvol*meanmass[i]/(-rhoc*Omegam*logsigmaM))
        if options.fType == 2: file.close()
	print "time takes= ",time.time()-tin

if options.fType == 3:
        print "\n####################### c-M relation #################################\n"
        def fitz(z, var):
          if var ==0:
            if z<=1: y= 1.08*(1+z)**0.27
            if z>1: y= 1.3*z**-0.19
          if var>0 and z<=1 and z>0: y=1.08*(1+z)**0.35 
          if var>0 and z==0: y=0.95
          return y
            # z0=1.08, z1=1.3, z2=1.14
        tin= time.time()      
        sodprofile= "sodproprtybins".rstrip()
        #halofileroot=fileroot+'/'+halofile
	cts_col= readheaders(options.cMFile, sodprofile,'sod_halo_bin_count')
	radius_col= readheaders(options.cMFile, sodprofile,'sod_halo_bin_radius')
        overden_col= readheaders(options.cMFile, sodprofile,'sod_halo_bin_rho_ratio')        
        print options.cMFile
        print overden_col
        print radius_col
        print cts_col
        print filenum
        print "\ndoing c-M calculation now.."
        conc, concerr, norm, massSO= conc_each_halo_lessmem(SOmass, SOradius, options.cMFile, sodprofile, filenum, SObinstart, cts_col, radius_col, overden_col)
        cmean, cmeanerr, meanmass, count, variance= get_bin_cmean(massSO, conc, concerr, bin_edges, massbins)

	file= open(options.outfile,'w')
        print>>file, "# SO Mass Msun/h\t# clusters\tcmean\tfrac err\tvariance\n"
        for i in range(0,massbins-1):  
            if count[i] > 3 and variance[i] !=0.0 and cmean[i] !=0:  print >>file, '%le\t%5d\t%3.4lf\t%3.4lf\t%3.4lf' %(meanmass[i], count[i],cmean[i]*fitz(zseek, model_no), (cmeanerr[i]**2+cmean[i]**2/count[i])**0.5, (variance[i]/cmean[i]**2-1.0)**0.5)
        #for i in np.arange (0,2.2,0.2): print i, fitz(i)*1.0/(1+i)**(0.71*0.89)
        file.close()
	print "time takes= ",time.time()-tin


"""
if case == 5:
     print "\n ######################### create mysql db########################\n"
     createcMdb('conc',massSO, conc, zseek,'cMdata')
  
  """ 

filelist = [ f for f in os.listdir("/project/projectdirs/hacc/PDACS/working/") if f.startswith("dataset") ]
for f in filelist:
    os.remove("/project/projectdirs/hacc/PDACS/working/" + f)
