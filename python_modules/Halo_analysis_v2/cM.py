from  numpy import *
from scipy  import *
from scipy.optimize import curve_fit

profnum=20
rhoc= 2.77536627e11; #Msun.h^2/Mpc^3
rmin=0.1 # minimum radius of the halo profile for the c-M fit
rmax= 1.8 # max radius

def lognfw (x,N,c):
   if c > 1.0 and c < 15.0 and N > 0.5 and N < 2:
     #y= log10(N*200.0/3.0) +3*log10(c) - log10(c * x)- 2*log10(1 +c*x) -log10(log(1+c)-c/(1+c)) 
     #y= N*200/3.0*c**3/(c*x)/(1+c*x)**2/(log(1+c)-c/(1+c))
     y= N*(log(1+c*x)-c*x/(1+c*x))/(log(1+c)-c/(1+c)) #/(4./3.*pi*x**3)#/rhoc
   else:
     y=1e10
   return y

def nfwdiff (x, delta, N,c):
   if c > 1.0 and c < 15.0 and N > 0.5 and N < 2:
     #y= log10(N*200.0/3.0) +3*log10(c) - log10(c * x)- 2*log10(1 +c*x) -log10(log(1+c)-c/(1+c)) 
     #y= N*200/3.0*c**3/(c*x)/(1+c*x)**2/(log(1+c)-c/(1+c))
     y= - N*(log(1+c*x)-c*x/(1+c*x))/(log(1+c)-c/(1+c)) #/(4./3.*pi*x**3)#/rhoc
     + N*(log(1+c*x*(1-delta))-c*x*(1-delta)/(1+c*x*(1-delta)))/(log(1+c)-c/(1+c))
   else:
     y=1e10
   return y

def get_bin_cmean(a, a2, a2_err, bin_edges, massbins):
   mean_a2= []
   mean_a2err= []
   mean_a=[]
   count=[]
   var_a=[]
   for i in range (0, massbins-1):
      b_start= bin_edges[i]
      b_end= bin_edges[i+1]
      ind_upper = nonzero(a >= b_start)[0]
      a_upper = a[ind_upper]
      a2_upper = a2[ind_upper]
      a2err_upper = a2_err[ind_upper]
      a_range = a_upper[nonzero(a_upper < b_end)[0]]
      a2_range = a2_upper[nonzero(a_upper < b_end)[0]]
      a2err_range = a2err_upper[nonzero(a_upper < b_end)[0]]
      if(len(a_range)!=0):
         meana= dot(a_range, a2_range)/sum(a_range)
         meanaerr = dot(a_range, a2err_range)/sum(a_range)
         vara= dot(a_range, a2_range**2)/sum(a_range)
         meanarange=mean(a_range)
      else: 
         meana, meanaerr,meanarange, vara=0., 0., 0., 0.
      mean_a2.append(meana)
      mean_a2err.append(meanaerr)
      mean_a.append(meanarange)
      var_a.append(vara)
      count.append(len(a_range))
   return array(mean_a2), array(mean_a2err), array(mean_a), array(count), array(var_a)

def conc_each_halo_lessmem(mass, r200, halofileroot, sodprofile,filenum, minmassSO, cts_col, radius_col, overden_col):
   conc, concerr, norm, MSO=[], [], [], []
   counter=0
   iter=0
   for i in range (0,filenum):
      infile= halofileroot 
#+"."+sodprofile+"."+str(i)
      radius, overden, cts= [], [], []    
      with open(infile,"r") as f:
        for line in f: 
          if not line.lstrip().startswith('#'):
             data = line.split()
             #print data,data[cts_col], data[radius_col], data[overden_col]
             #exit()
             if float(data[radius_col])/r200[iter] > rmin and float(data[radius_col])/r200[iter]  < rmax:
                radius.append(float(data[radius_col]))
                overden.append(float(data[overden_col]))
                cts.append(int(data[cts_col]))
             counter= counter+1
             remainder,q=modf(float(counter)/float(profnum))
             if remainder==0:
               if mass[iter]> minmassSO:
                  Delta= array(overden)*rhoc*(4./3.*pi*(array(radius))**3)/mass[iter]
                  r=radius/r200[iter]
                
                  
                  #yerr=cumsum(cts[0:profnum])**(-0.5)
                  #print radius/r200[iter]
                  #for i=0
                  #if radius/r200[iter] > 0.1:
                  
                  #   print r
                  #exit()                          
                  #print r, Delta
                  #if iter>2: exit()
                  N0, c0=1.0, 5.0
                  pfinal, covar=curve_fit(lognfw, r, Delta, p0=[N0, c0], maxfev=10000)
                  conc.append(pfinal[1])
                  concerr.append(covar[1][1]**0.5)
                  norm.append(pfinal[0]) 
                  MSO.append(mass[iter])
               radius, overden, cts= [], [], []
               iter=iter+1 
               #r1,q1=modf(float(iter)/100)
               #if r1==0 : print iter,' fitting halos done'               
   return array(conc), array(concerr), array(norm), array(MSO)

def conc_each_halo_lessmem_diff(mass, r200, halofileroot, sodprofile,filenum, minmassSO, cts_col, radius_col, overden_col):
   conc, concerr, norm, MSO=[], [], [], []
   counter=0
   iter=0
   for i in range (0,filenum):
      infile= halofileroot
#+"."+sodprofile+"."+str(i)
      radius, overden, cts= [], [], []    
      with open(infile,"r") as f:
        for line in f: 
          if not line.lstrip().startswith('#'):
             data = line.split()
             #print data[cts_col], data[radius_col], data[overden_col-2]
             #exit()
             if float(data[radius_col])/r200[iter] > rmin and float(data[radius_col])/r200[iter]  < rmax:
                radius.append(float(data[radius_col]))
                overden.append(float(data[overden_col-2]))
                cts.append(int(data[cts_col]))
             counter= counter+1
             remainder,q=modf(float(counter)/float(profnum))
             if remainder==0:
               if mass[iter]> minmassSO:
                  Delta= array(overden)*rhoc*(4./3.*pi*(array(radius))**3)/mass[iter]
                  r=radius/r200[iter]
                  deltar= radius[9]/radius[10]
                  #print deltar
                  #if counter==100: exit()
                  #yerr=cumsum(cts[0:profnum])**(-0.5)
                  #print radius/r200[iter]
                  #for i=0
                  #if radius/r200[iter] > 0.1:
                  
                  #   print r
                  #exit()                          
                  #print r, Delta
                  #if iter>2: exit()
                  N0, c0=1.0, 5.0
                  pfinal, covar=curve_fit(lambda x, N, c: nfwdiff(x, deltar, N,c), r, Delta, p0=[N0, c0], maxfev=10000)
                  conc.append(pfinal[1])
                  concerr.append(covar[1][1]**0.5)
                  norm.append(pfinal[0]) 
                  MSO.append(mass[iter])
               radius, overden, cts= [], [], []
               iter=iter+1 
               #r1,q1=modf(float(iter)/100)
               #if r1==0 : print iter,' fitting halos done'
               
   return array(conc), array(concerr), array(norm), array(MSO)


def readhaloprof(halofileroot, filenum, profnum):
   radius= []
   overden= []
   cts=[]
   for i in range (0,filenum):
      infile= halofileroot
#+".sodprofile."+str(i)
      for line in open(infile,'r'):
           data = line.split(' ')
           #print data
           radius.append(float(data[radius_col]))
           overden.append(float(data[overden_col]))
           cts.append(int(data[cts_col]))
                
   return array(radius), array(overden), array(cts)

def conc_each_halo(mass, r200, radius, overden, cts, minmassSO):
   conc=[]
   concerr=[]
   norm=[]
   MSO=[]
   for i in range (0, len(mass)):
     if(mass[i]> minmassSO):
        r=(radius[i*profnum:(i+1)*profnum]) 
        Delta=(overden[i*profnum:(i+1)*profnum])*rhoc*(4./3.*pi*r**3)/mass[i]
        csum = cumsum(cts[i*profnum:(i+1)*profnum])**0.5
        #if csum[0] !=0: yerr=csum**0.5
        #else:           yerr= csum+1
        #else: yerr.all=1.0
        #print yerr, 

        r=r/r200[i]   
        #print r,Delta
        #if i>2: exit()
        N0, c0=1.0, 5.0
        pfinal, covar=curve_fit(lognfw, r, Delta, p0=[N0, c0],maxfev=10000)
        conc.append(pfinal[1])
        concerr.append(covar[1][1]**0.5)
        norm.append(pfinal[0])
        MSO.append(mass[i])
        # print i
        #exit() 
   return array(conc), array(concerr), array(norm), array(MSO)

