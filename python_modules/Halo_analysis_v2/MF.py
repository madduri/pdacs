from numpy import *
from sys import *
from scipy import *

def read_input_params (infile):
   inputparams= []
   with open(infile,"r") as f:
     for line in f:
       if not line.lstrip().startswith('#'):
          inputparams.append(line)
   return array(inputparams)
  
def readheaders(halofileroot,  properties, coltype):
     #infile= halofileroot+"."+properties+".0"
     infile = halofileroot
     with open(infile,"r") as f:
       for line in f:
         if line.lstrip().startswith('#') and (line.find('halo_count') > 0 or line.find('sod_halo_bin')>0):
            hdr= line.split()
            return hdr.index(coltype)-1 
     
def readhalofiles(halofileroot, properties, count_col, filenum):
   count= []
   for i in range (0,filenum):
      #infile= halofileroot+"."+properties+"."+str(i)
      infile=halofileroot
      with open(infile,"r") as f:
        for line in f:
           if not line.lstrip().startswith('#'):
                 data = line.split()
                 #print data,data[count_col]
                 #exit()
                 count.append(int(data[count_col]))              
   return array(count)

def get_bin_mean(a,bin_edges, massbins):
   mean_mass= []
   for i in range (0, massbins-1):
      b_start= bin_edges[i]
      b_end= bin_edges[i+1]
      ind_upper = nonzero(a >= b_start)[0]
      a_upper = a[ind_upper]
      a_range = a_upper[nonzero(a_upper < b_end)[0]]
      if(len(a_range)!=0): meana= mean(a_range)
      else: meana=0.5*(b_start+b_end)
      mean_mass.append(meana)
   return array(mean_mass)
     

def logrange(start, end, binnum):
   b_edge=[]
   binsize=log(end/start)/(binnum-1)
   for i in range(0,binnum):
      b_edge.append(start*exp(i*binsize))     
   return array(b_edge), binsize

def MF_fit(x, z):
  A=0.333/(1+z)**0.11
  a1= 0.788/(1+z)**0.01
  p= 0.807
  q= 1.795
  if z==0: delc=1.674
  elif z==1: delc=1.684
  elif z==2: delc= 1.686
  else: delc= 1.684

  return A*(2/pi)**0.5*exp(-a1*delc**2/2/x**2)*(1+(x**2/a1/delc**2)**p)*(delc/x*a1**0.5)**q

def calc_mf(b_start, b_end, massbins, mass):
    bins, binsize= logrange(b_start, b_end, massbins) 
    hist,bin_edges=histogram(mass,bins=bins)
    meanmass= get_bin_mean(mass, bin_edges, massbins)
    return meanmass, hist, binsize, bin_edges
