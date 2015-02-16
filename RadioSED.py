#!/usr/bin/env python


import os,sys

#import pyfits
#from numpy import *
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import patches
from matplotlib.font_manager import FontProperties
from pylab import title,xlabel,ylabel,text

import Cosmology

# usage:
#                redshift
#    RadioSED.py 1.25

# constants
c = 2.9979E10       # cm/s
mu = 1E-4           # cm
Mpc = 3.08568025E24 # cm

Jy = 1.0E-23        # erg/s/cm^2/Hz
Lsun = 3.839E33     # erg/s

z = float(sys.argv[1])

# ALMA
# 23.0278 deg S, 67.7548 deg W


# Table 1
# https://almascience.nrao.edu/proposing/call-for-proposals/capabilities
# updated 11-11-13

# Notes for Table 1:
#   These are the nominal frequency ranges for continuum observations.
#   Observations of spectral lines that are within about 0.2 GHz of a band edge
#   are not possible at present in Frequency Division Mode (FDM, see section
#   A.6.1), because of the responses of the spectral edge filters implemented
#   in the correlators.

# 1 hr sensitivity (rms) for 0.1" resolution for 34 12m Antennas at -55 dec
# 84 - 116 GHz
# 0.01058 - 0.03470 mJy
# 1 GHz sampling
sens3 = [0.01058, 0.01329, 0.01326, 0.01323, 0.01323, 0.01323, 0.01324,
         0.01326, 0.01328, 0.01332, 0.01349, 0.01340, 0.01350, 0.01351,
         0.01356, 0.01363, 0.01371, 0.01381, 0.01397, 0.01257, 0.01271,
         0.01276, 0.01290, 0.01306, 0.01328, 0.01355, 0.01394, 0.01465,
         0.01518, 0.01637, 0.02038, 0.02474, 0.03470]

# 125 - 163 GHz
# 0.01760 - 0.01529 mJy
# 5 GHz sampling
sens4 = [0.01760, 0.01521, 0.01488, 0.01496, 0.01549, 0.01419, 0.01448,
         0.01491, 0.01529]

# 211 - 275 GHz
# 0.01861 - 0.01791 mJy
# 5 GHz sampling
sens6 = [0.01861, 0.01773, 0.01664, 0.01686, 0.01707, 0.01769, 0.01780,
         0.01768, 0.02280, 0.01824, 0.01850, 0.01928, 0.01752, 0.01791]

# 275 - 373 GHz
# 0.01791 - 0.06679 mJy
# 5 GHz sampling
sens7 = [0.01791, 0.02428, 0.02464, 0.02296, 0.02338, 0.02410, 0.02471,
         0.02625, 0.02799, 0.04402, 0.21638, 0.03428, 0.02745, 0.02649,
         0.02682, 0.02777, 0.03432, 0.03269, 0.04140, 0.07639, 0.06679]

# 385 - 500 GHz
# 0.19437 - 0.12351 mJy
# 5 GHz sampling
sens8 = [0.19437, 0.07378, 0.06761, 0.06566, 0.05885, 0.05998, 0.06559,
         0.09312, 55080210.14793, 0.10294, 0.08177, 0.21751, 0.81174,
         5.14005, 0.12489, 0.08427, 0.08100, 0.14609, 0.58080, 0.10704,
         0.18370, 0.18150, 0.11436, 0.12351]

# 420 - 430 GHz
# 1 GHz sampling
# 0.09312, 0.12351, 0.18918, 0.67878, 88.33793, 55080210.14793, 3.34911,
# 0.32239, 0.16061, 0.12086, 0.10294

# 445-455 GHz
# 1 GHz sampling
# 0.81174, 4.93033, 384.71663, 155941.38758, 448.50170, 5.14005, 0.81965,
# 0.33891, 0.20648, 0.15247, 0.12489

#602    , 720
#0.62789, 0.94408
sens9 = [0.62789, 0.54369, 0.51763, 0.56534, 181.37324, 0.63486, 0.33324,
         0.28805, 0.27535, 0.27329, 0.24727, 0.84885, 0.25983, 0.24468,
         0.27923, 0.24680, 0.25416, 0.26624, 0.28386, 0.31728, 0.34181,
         0.38903, 0.48843, 21.81105, 0.94408]

# frequency in GHz
alma = {3:{"nu0":84 ,"nu1":116,"sens":sens3},
        4:{"nu0":125,"nu1":163,"sens":sens4},
        6:{"nu0":211,"nu1":275,"sens":sens6},
        7:{"nu0":275,"nu1":373,"sens":sens7},
        8:{"nu0":385,"nu1":500,"sens":sens8},
        9:{"nu0":602,"nu1":720,"sens":sens9}}

vla  = {"L band" :{"nu0":1 ,  "nu1":2},     # low frequencies
        "S band" :{"nu0":2 ,  "nu1":4},
        "C band" :{"nu0":4 ,  "nu1":8},
        "X band" :{"nu0":8 ,  "nu1":12},
        "Ku band":{"nu0":12,  "nu1":18},    # high frequencies
        "K band" :{"nu0":18,  "nu1":26.5},
        "Ka band":{"nu0":26.5,"nu1":40},
        "Q band" :{"nu0":40,  "nu1":50}}

# v = 0     GHz
CO = {"1-0":115.27120,
      "2-1":230.53800,
      "3-2":345.79599,
      "4-3":461.04077,
      "5-4":576.26793,
      "6-5":691.47308,
      "7-6":806.65180,
      "8-7":921.79970,}

# variables
wmin = 10
#wmax = 4E3  # ALMA
wmax = 5E5 # ALMA+VLA

write = 0
if "--write" in sys.argv: write = 1

# Calculate z at 10 Mpc for templates
c1 = Cosmology.Cosmology()
Dlt = 10            # Mpc
zt = c1.Dl2z(Dlt)

# Calculate Dl for galaxy at z
Dl = c1.Dl(z)            # Mpc
print Dl

tab = "sed_templates/table4_avgtemplates.txt"
table = np.loadtxt(tab,skiprows=26)      # table4_avgtemplates.txt

wavT = table[:,0]          # um
fluxT = table[:,1:]        # Jy

fmin = 1000
fmax = -1000
mir_min = 1000
mir_max = -1000

# new wmin
wmin = np.min((1+z)*wavT)
print wavT[0]

# avgtemplates
templates = [9.75, 10.00, 10.25, 10.50, 10.75, 11.00, 11.25, 11.50, 11.75, 12.00, 12.25, 12.50, 12.75, 13.00]
#L = SFR/4.5E-44            # erg/s   FIR  Kennicutt Review

SFR = 4.5E-44*10**(np.array(templates))*Lsun
print templates
print SFR

# localseds
redshifts = [0.015938,0.010807,0.009354,0.008342,0.073317,0.083000,0.0429,0.077760,0.018126,0.010781,0.012999]

#fig = plt.figure(figsize=(7,6))
fig = plt.figure(figsize=(10,6))
p = fig.add_subplot(111)
#p.set_axis_bgcolor("k")

#colormap = plt.cm.gist_ncar
num_plots = len(templates)
colormap = plt.cm.Greys
colors = [colormap(i) for i in np.linspace(0.5, 0.9, num_plots)]

# Galaxy SEDs
for j,galaxy in enumerate(np.transpose(fluxT)):
    Fnu_l = galaxy*Jy                 # erg/s/cm^2/Hz
    # calculate Lnu for the template at 10 Mpc
    Lnu = 4*np.pi*(Dlt*Mpc)**2*Fnu_l/((1+zt)) 
    # calculate Fnu of the temaplate at z
    Fnu = Lnu/(4*np.pi*(Dl*Mpc)**2)*(1+z)/Jy   # Jy


    if min(Fnu) < fmin and min(Fnu) != 0.0:
        fmin = min(Fnu)

    if max(Fnu) > fmax:
        fmax = max(Fnu)

    #filt = np.greater_equal(wavT,wmin)*np.less_equal(wavT,10)
    filt = np.greater_equal(wavT,wmin)*np.less_equal(wavT,wmax)
    mir = np.compress(filt,Fnu,0)

    if min(mir) < mir_min and min(mir) != 0.0:
        mir_min = min(mir)

    if max(mir) > mir_max:
        mir_max = max(mir)

    #p.plot(wavT*(1+z),Fnu)
    p.plot(wavT*(1+z),Fnu,c=colors[j])

print fmin,fmax
print fmin/100.,fmax*10

trans = CO.keys()
trans.sort()
print trans

Nu0 = np.array([alma[b]["nu0"] for b in alma]+[vla[b]["nu0"] for b in vla])
Nu1 = np.array([alma[b]["nu1"] for b in alma]+[vla[b]["nu1"] for b in vla])

# CO lines
for t in trans:
    nu_emit   = CO[t]*1E9   # Hz
    nu_obs    = nu_emit/(1+z)
    lamb_obs  = c/nu_obs    
    lamb_emit = c/nu_emit

    # determine if line is in ALMA band
    #   count the number of positive signs
    #     equal   == outside of the ALMA band
    #     unequal == inside the ALAM band

    #         band edge --> not accessible for lines
    Vec0 = Nu0 + 0.2 - nu_obs/1E9  # minimum frequency for an ALMA band
    Vec1 = Nu1 - 0.2 - nu_obs/1E9  # maximum frequency ...

    n0 = sum(np.greater(Vec0,0))
    n1 = sum(np.greater(Vec1,0))
  
    if n0 == n1: color = "r"
    else:        color = "g"

    p.plot([lamb_obs/mu,lamb_obs/mu],[fmin/100.,fmax*10],"--",c=color,lw=2)

    p.text(0.95*lamb_obs/mu,fmax/10.,"CO (%s)" % t,fontsize=8,
        rotation='vertical',
        horizontalalignment='center',
        verticalalignment='bottom',)

# ALMA bands
for b in alma:
    nu0 = alma[b]["nu0"]
    nu1 = alma[b]["nu1"]
    sens = np.array(alma[b]["sens"])*1e-3*5

    w1 = c/(nu0*1E9)/mu
    w0 = c/(nu1*1E9)/mu

    xy = (w0, fmin/100.)
    rect = patches.Rectangle(xy,w1-w0,fmax*10,alpha=0.2)
    p.add_artist(rect)

    p.text((w0+w1)/2,fmax,"Band %i" % b,fontsize=8,
        rotation='vertical',
        horizontalalignment='center',
        verticalalignment='bottom',)

    p.plot(np.linspace(w0,w1,len(sens)),sens,"--",c="k",lw=2) # ~1hr

for b in vla:
    nu0 = vla[b]["nu0"]
    nu1 = vla[b]["nu1"]
    #sens = alma[b]["sens"]

    w1 = c/(nu0*1E9)/mu
    w0 = c/(nu1*1E9)/mu

    xy = (w0, fmin/100.)
    rect = patches.Rectangle(xy,w1-w0,fmax*10,alpha=0.2)
    p.add_artist(rect)


    if (w0+w1)/2 > wmin and (w0+w1)/2 < wmax:
        p.text((w0+w1)/2,fmax,b,fontsize=8,
            rotation='vertical',
            horizontalalignment='center',
            verticalalignment='bottom',)

p.plot([wmin,wmax],[10E-3,10E-3],"-.",c="k") # 10 mJy source
#p.plot([wmin,wmax],[1E-3,1E-3],"--",c="k") # 1 mJy source
p.plot([wmin,wmax],[1E-4,1E-4],"--",c="k") # 0.1 mJy source  ~2hr


#p.text(5.*(1+z),mir_min/6.,"$10^{%.1f}$ L$_{\odot}$" % templates[0],fontsize=8)
#p.text(5.*(1+z),mir_max,"$10^{%.1f}$ L$_{\odot}$" % templates[-1],fontsize=8)
p.text(0.1,0.25,"$10^{%.1f}$ L$_{\odot}$" % templates[0],fontsize=12,transform=p.transAxes)
p.text(0.1,0.8,"$10^{%.1f}$ L$_{\odot}$" % templates[-1],fontsize=12,transform=p.transAxes)


ax = p.twiny()
ax.set_xscale("log")
ax.set_xlim(c/(wmin*mu)/1E9,c/(wmax*mu)/1E9) # in GHz
ax.set_xlabel("Frequency (GHz)")

p.set_xlim(wmin,wmax)
p.set_ylim(fmin/100.,fmax*10)
p.set_xscale("log")
p.set_yscale("log")
p.set_xlabel("Wavelength ($\mu$m)")
p.set_ylabel("Flux (Jy)")

p.text(0.05,0.95,"z = %.2f" % z,fontsize=12,transform=p.transAxes)

if write:
    #fig.savefig("RadioSED.eps")
    fig.savefig("RadioSED.pdf")
else:
    plt.show()
