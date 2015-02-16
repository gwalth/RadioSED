
#from math import sqrt
#from Numeric import *
#import pgplot
#from scipy import *

from math import sqrt
from numpy import *
from scipy import *
from scipy import interpolate,integrate

# constants
c = 2.9979E5             # km/s        9.46E24 m/Gyr

# preset values
h = 0.7
Om_r = 0.0  # Om_r = 1E-5
Om_m = 0.3
Om_l = 0.7
Om_K = 0.0
z1 = 0.0
z2 = 13.0
zint = 0.001

class Cosmology:
    def __init__(self, h=h, Om_r=Om_r, Om_K=Om_K, Om_m=Om_m, Om_l=Om_l,
                 z1=z1, z2=z2, zint=zint):
        self.h = h
        self.Ho = 100*h               # km/s/Mpc 
        self.Dh = c/self.Ho
        
        self.Om_r = Om_r
        self.Om_m = Om_m
        self.Om_l = Om_l
        # latest change, not sure if this is bad
        if Om_K:
           self.Om_K = Om_K
        else:
           self.Om_K = 1.0 - Om_r - Om_m - Om_l

        self.z = arange(z1,z2,zint)
        self.z1 = z1
        self.z1 = z2
        self.zint = zint
        # cheat, creating lookup table to call function only once
        self.rz0 = self.rz()

    def Sr(self,r):       
        Rc = None
        if self.Om_K == 0:      # K = 0            # Flat
            S = r
        else:
            if self.Om_K > 0:   # K < 0            # Hyperbolic
                Rc = sqrt(c**2/(self.Ho**2*self.Om_K))  # Mpc
                S = Rc*sinh(r/Rc)
            elif Om_K < 0:      # K > 0            # Spherical
                self.Rc = sqrt(-c**2/(self.Ho**2*self.Om_K)) # Mpc
                S = Rc*sin(r/Rc)
        return S

    def rz(self):
        # comoving distance
        # the loop is slow, need a faster way to compute
        return [integrate.quad(self.rdz,0,z1)[0] for z1 in self.z]
        # maybe a way to avoid the loop, use trapz
        #return integrate.trapz(c/self.Hz(self.z))

    def rdz(self,z):
        return c/self.Hz(z)
    
    def dz(self,z):
        return 1/((1+z)*self.Hz(z))

    def Hz(self,z):
        return self.Ho*sqrt(self.Om_r*(1+z)**4 + self.Om_m*(1+z)**3 \
                            + self.Om_K*(1+z)**2 + self.Om_l)

    def Da_arr(self):                    
        #return self.Sr(self.rz())/(1+self.z)
        return self.Sr(self.rz0)/(1+self.z)

    def Dl_arr(self):
        #return self.Sr(self.rz())*(1+self.z)
        return self.Sr(self.rz0)*(1+self.z)
    
    def Da(self,z0): 
        # angular distance array
        Da_ = interpolate.interp1d(self.z,self.Da_arr())
        return Da_(z0)

    def Dl(self,z0): 
        # luminosity distance array
        Dl_ = interpolate.interp1d(self.z,self.Dl_arr())
        return Dl_(z0)

    def Dl2z(self,Dl):
        z_ = interpolate.interp1d(self.Dl_arr(),self.z)
        return z_(Dl)

    def Da2z(self,Da):
        z_ = interpolate.interp1d(self.Da_arr(),self.z)
        return z_(Da)

    def ti(self,z1,z2):
        # lookback  0 to z
        # age       z to inf
        return integrate.quad(self.dz,z1,z2)[0]
        

def test():
    c = Cosmology()
    print c.Da(0.6107)
    print c.Dl(0.6107)



