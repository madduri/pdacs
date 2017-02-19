from scipy.integrate import odeint 
#from scipy.interpolate import interp1d
#from pylab import *
import numpy as np

"""
solve the cosmological growth factor ODE:
D'' +N1(a)*D +N2(a)*D'=0
where N1(a)= 3/2*a^-2* X(a)/(1+X(a))
      N2(a)= -3/(2a)*(1-w(a)/(1+X(a)))
      X(a)= Om*1/(E^2(a)*a^3-Om)
      E^2(a)=H/H0=Om*a^-3+(1-Om)exp(-3*int_1^a' d ln a (1+wt) )
"""

def deriv(y,a, Om, Ol, wt):
 # return derivatives of the array y 
  return np.array([ y[1], N1(a, Om, Ol, wt)*y[0]+N2(a, Om, Ol, wt)*y[1]])

def N1(a, Om, Ol, wt):
  # return the func N1
  return  3./2./a**2*X(a, Om, Ol, wt)/(1+X(a, Om, Ol, wt))

def N2(a, Om, Ol, wt):
  #return the func N2
  return -3./2./a*(1-wt/(1+X(a, Om, Ol, wt)))

def X(a, Om, Ol, wt):
 # return the func X
 return Om/Ol*a**(3*wt)

def Da_ode(Om, Ol, wt):
  ainitial=1e-4
  afinal=1.0
  size=10000
  scalefac = np.linspace(ainitial,1.0,size) 
  yinit = np.array([ainitial, 0.0]) # initial values y = 
  yf=odeint(deriv,yinit,scalefac, args=(Om, Ol, wt))
  D = yf[:,0]/yf[size-1,0]
  logD = yf[:,1]*scalefac/yf[:,0]
  return (scalefac, D, logD)

def interp_D(zseek, ascale, D0, D1):
  aseek= 1/(1+zseek)
  Dseek=np.interp(aseek, ascale, D0)
  logDseek= np.interp(aseek, ascale, D1)
  return (Dseek, logDseek)

#def plot_growth(scalefac, D, logD):
#  figure() 
#  plot(scalefac, D) # y[:,0] is the first column of y 
#  plot(scalefac, logD)
#  xlabel("a") 
#  ylabel("D(a)/D(1); log D/ log a") 
#  show()

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 5:
        print "usage: python growth_ode.py Omega_m Omega_l w0 redshift"
        exit()

    scalefac, D0, D1= Da_ode(float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]))
    zseek= float(sys.argv[4])
    Ds, lDs= interp_D(zseek,scalefac,D0, D1)
    print "redshift D(z) dlog D/d log a (z)"
    print zseek, Ds, lDs
    #plot_growth(scalefac, D0, D1)
