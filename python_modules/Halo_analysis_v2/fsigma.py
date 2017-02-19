from numpy import *
from sys import *
from scipy import *
from scipy.integrate import *

rhoc= 2.77536627e11; #Msun.h^2/Mpc^3

def read_Tk(fileTk):
   print fileTk
   k=[]
   Tk=[]
   for line in open(fileTk,'r'):
     data = line.split()
     k.append(float(data[0]))
     Tk.append(float(data[6])) #6
   k=array(k)
   Tk=array(Tk)
   return k,Tk

def Pk_norm(sigma8, ns, k, Tk, hubble):
  R=8#/hubble
  y= k**(ns+2)*Tk**2*(3*(sin(k*R)-k*R*cos(k*R))/(k*R)**3)**2
  result= simps(y,k,even='avg')
  return sigma8**2*2*pi**2/result   
     
def sigmam(Om, ns, M, N, k, Tk, hubble):
  R= (3*M/(4*pi*Om*rhoc))**(1./3.)#/hubble
  y= N/(2*pi**2)*k**(ns+2)*Tk**2*(3*(sin(k*R)-k*R*cos(k*R))/(k*R)**3)**2
  return simps(y,k,even='avg')**0.5

def logsigm(Om, ns, M, N, k, Tk, sigmam, hubble):
   R= (3*M/(4*pi*Om*rhoc))**(1./3.)#/hubble
   #y= N*k**(ns+1)*Tk**2*3*(sin(k*R)-k*R*cos(k*R))/(k*R)**3*(sin(k*R)+3*cos(k*R)/(k*R)-3*sin(k*R)/(k*R)**2)
   WkR= 3*(sin(k*R)-k*R*cos(k*R))/(k*R)**3
   y=N*k**(ns-1)*Tk**2*WkR*(((k*R)**2-3)*sin(k*R)+ 3*k*R*cos(k*R))
   return 1/(2*pi**2*R**3)*simps(y,k,even='avg')/sigmam**2
   
if __name__ == "__main__":
    import sys
    infile= sys.argv[1]
    ns=0.97
    sigma8=0.8
    h=0.72
    Om=0.25 
    M=1e15
    k,Tk= read_Tk(infile)
    N= Pk_norm(sigma8, ns, k, Tk, h)

    sigmaM= sigmam(Om, ns, M, N, k, Tk, h)
    logsigmaM= logsigm(Om, ns, M, N, k, Tk, sigmaM, h)
    print M, 1/sigmaM, logsigmaM
