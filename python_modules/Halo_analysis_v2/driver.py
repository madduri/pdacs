import sys
from MF import *
from cM import *
from fsigma import *
from growth_ode import *
#from createcMdb import *

import time

rhoc= 2.77536627e11; #Msun.h^2/Mpc^3

tin= time.time()
infile= sys.argv[1]
inputparams= read_input_params(infile) 

case= int(inputparams[0])
massbins= int(inputparams[1])
filenum= int(inputparams[2])
zseek= float(inputparams[3])
if int(inputparams[4])==1:
  model_no=int(inputparams[5])
  cosmoparams= read_input_params("model_param.txt")  
  tmp=cosmoparams[model_no].split()
  hubble= float(tmp[8])
  Omegam=float(tmp[6])
  ns= float(tmp[3])
  sigma8= float(tmp[4])
  w0= float(tmp[5])
else:
	model_no=0
	hubble= float(inputparams[6])
	Omegam= float(inputparams[7])
	ns= float(inputparams[8])
	sigma8= float(inputparams[9])
	w0= float(inputparams[10])
print hubble, Omegam, ns, sigma8, w0
#exit()
boxsize= float(inputparams[11])
particlenum= int(inputparams[12])
minparticle= int(inputparams[13])
Delta= float(inputparams[14])
fileroot= inputparams[15].rstrip()
halofile= inputparams[16].rstrip()
Tkfile= inputparams[17].rstrip()

massres= Omegam*rhoc*(boxsize/particlenum)**3 #Msun/h
simvol= (boxsize)**3 #Mpc^3/h^3
print 'mass res= %le' % (massres), 'Msun/h'    
    
k,Tk= read_Tk(Tkfile)
N= Pk_norm(sigma8, ns, k, Tk, hubble)
scalefac, D0, D1= Da_ode(Omegam, 1-Omegam, w0)
Ds, lDs= interp_D(zseek,scalefac,D0, D1)
print "D(z) and log D/log z(z)= ", Ds, lDs, 'at z= ',zseek
print "time takes= ",time.time()-tin

if case == 0 or case == 2 or case==4:
        print "\n################## FOF MF ##########################################\n"
        fofproperties= inputparams[18].rstrip()
        halofileroot=fileroot+'/'+halofile
        fofcount_col= readheaders(halofileroot, fofproperties,'fof_halo_count')
        FOFcount= readhalofiles(halofileroot, fofproperties, fofcount_col,filenum)
        
        print "\ntotal # of fof clusters= ",len(FOFcount)
        FOFbinstart=min(FOFcount)*massres
        FOFbinend=max(FOFcount)*massres
        FOFmass=FOFcount*(1.0-1.0/FOFcount**0.65)*massres
        print "minFOF maxFOF"
        print "# particles= ",min(FOFcount), max(FOFcount)
	print "mass Msun/h= ", '%le %le' %(FOFbinstart, FOFbinend)
        file= open(halofile+'.fof_mf','w')
        print >>file, "# FOF Mass Msun/h\t# clusters\tdn/dln M (h/Mpc)^3\tfrac err\t1/sigma(M)\tf(sigma)\tfsigma_fit\t\n"
        meanmass, hist, binsize, bin_edges= calc_mf(FOFbinstart, FOFbinend, massbins, FOFmass)
        for i in range(0,massbins-1):  
           sigmaM= sigmam(Omegam, ns, meanmass[i], N, k, Tk, hubble)
           logsigmaM= logsigm(Omegam, ns, meanmass[i], N, k, Tk, sigmaM, hubble)
           sigmaM= Ds*sigmaM
           if hist[i] >0:  print >> file, '%le\t%5d\t%7.4le\t%7.4lf\t%7.4lf\t%7.4lf\t%7.4lf' %(meanmass[i], hist[i], hist[i]/binsize/simvol,1/hist[i]**0.5, 1/sigmaM, hist[i]/binsize/simvol*meanmass[i]/(-rhoc*Omegam*logsigmaM), MF_fit(sigmaM, zseek))
	file.close()

if case >= 1 and case<=5:
        print "\n##################### SO MF ########################################\n"
        tin= time.time()
        sodproperties= inputparams[19].rstrip()
        halofileroot=fileroot+'/'+halofile
        halofileroot="../working/dataset_2637.dat"
        socount_col= readheaders(halofileroot, sodproperties,'sod_halo_count')
        SOcount= readhalofiles(halofileroot, sodproperties, socount_col, filenum)   
	print "\ntotal # of SO clusters= ",len(SOcount)
        SObinstart=max(min(SOcount), minparticle)*massres
        SObinend=max(SOcount)*massres
        SOmass=SOcount*massres
        SOradius= (3.0*SOmass/4.0/pi/rhoc/Delta)**(1.0/3.0)
	print "minSO maxSO"
        print "# particles= ",min(SOcount), max(SOcount)
	print "mass Msun/h= ", '%le %le' %(SObinstart, SObinend)
        file= open(halofile+'.sod_mf','w')
        print>> file, "# SO Mass Msun/h\t# clusters\tdn/dln M(h/Mpc)^3\tfrac err\t1/sigma(M)\tf(sigma)\n"
        meanmass, hist, binsize, bin_edges= calc_mf(SObinstart, SObinend, massbins, SOmass)
        for i in range(0,massbins-1):
           sigmaM= sigmam(Omegam, ns, meanmass[i], N, k, Tk, hubble)
           logsigmaM= logsigm(Omegam, ns, meanmass[i], N, k, Tk, sigmaM, hubble)
	   sigmaM= Ds*sigmaM
           if hist[i] >2:  print>> file, '%le\t%5d\t%7.4le\t%7.4lf\t%7.4lf\t%7.4lf' %(meanmass[i], hist[i], hist[i]/binsize/simvol,1/hist[i]**0.5, 1/sigmaM, hist[i]/binsize/simvol*meanmass[i]/(-rhoc*Omegam*logsigmaM))
        file.close()
	print "time takes= ",time.time()-tin

if case >= 3:
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
        sodprofile= inputparams[20].rstrip()
        halofileroot=fileroot+'/'+halofile
        halofileroot="../working/dataset_2640.dat"
	cts_col= readheaders(halofileroot, sodprofile,'sod_halo_bin_count')
	radius_col= readheaders(halofileroot, sodprofile,'sod_halo_bin_radius')
        overden_col= readheaders(halofileroot, sodprofile,'sod_halo_bin_rho_ratio')        
        print "\ndoing c-M calculation now.."
        print cts_col
        print radius_col
        print overden_col 
        conc, concerr, norm, massSO= conc_each_halo_lessmem(SOmass, SOradius, halofileroot, sodprofile, filenum, SObinstart, cts_col, radius_col, overden_col)
        cmean, cmeanerr, meanmass, count, variance= get_bin_cmean(massSO, conc, concerr, bin_edges, massbins)

	file= open(halofile+'.cM','w')
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


