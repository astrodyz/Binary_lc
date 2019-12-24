#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 10 14:22:47 2019

@author: dyz
"""

import phoebe
from phoebe import u
import numpy as np
import matplotlib.pyplot as plt

logger = phoebe.logger(clevel='WARNING')

kepler = phoebe.get_passband('Kepler:mean')
tess = phoebe.get_passband('TESS:default')

"""
calculate the two passband flux
"""
teffs = np.linspace(3500, 8000, 100)
fig=plt.figure()
plt.xlabel('Temperature [K]')
plt.ylabel('Inorm [W/m^2/A]')
plt.plot(teffs, kepler.Inorm(teffs, atm='blackbody', ld_func='linear', ld_coeffs=[0.0]), label='Kepler')
plt.plot(teffs, tess.Inorm(teffs, atm='blackbody', ld_func='linear', ld_coeffs=[0.0]), label='TESS')
plt.legend(loc='lower right')

plt.savefig("two_passband_flux.eps")

"""
get binary flux, default teff1=6000K, teff2=6000K,q=1,ecc=0, model=ck2004,atlas9
"""
def get_flux_value_kepler(teff):
    teffs = np.linspace(4000, 20000, 17)
    fluxes=kepler.Inorm(teffs, atm='blackbody', ld_func='linear', ld_coeffs=[0.0])
    for i in range(0,len(fluxes)):
        if teffs[i]==teff:
            return fluxes[i]
 
def get_flux_value_tess(teff):
    teffs = np.linspace(4000, 20000, 17)
    fluxes=tess.Inorm(teffs, atm='blackbody', ld_func='linear', ld_coeffs=[0.0])
    for i in range(0,len(fluxes)):
        if teffs[i]==teff:
            return fluxes[i]       
        
def get_binary_flux_kepler(teff1,teff2,q,ecc):
    base_flux=get_flux_value_kepler(teff1)
    b = phoebe.default_binary()
    b.set_value('q',value=q)
    b.set_value('ecc',value=ecc)
    b.set_value('teff', component='primary', value=teff1)
#    b.set_value('r',component='secondary',value=0.5)
    b.set_value('teff', component='secondary', value=teff2)
    b.set_value('q',value=0.75)#b.add_compute('lc',atm@primary='ck2004',atm@second='ck2004')
    b.add_dataset('lc',dataset='lc01',atm='blackbody',ld_func='interp',times = phoebe.linspace(0,2,101),passband='Kepler:mean',intens_weighting='energy')
#    b.add_dataset('lc',dataset='lc01',ld_func='interp',times = phoebe.linspace(0,1,101),passband='Kepler:mean',intens_weighting='energy')
    b.run_compute(irrad_method='none')
#    b['lc01@model'].plot(x='phases',c='blue')
#    b.savefig('lc_kepler.eps')
    times1=b['value@times@lc01@latest@model']
    fluxes1=b['value@fluxes@lc01@latest@model']*base_flux
    return times1,fluxes1


def get_binary_flux_tess(teff1,teff2,q,ecc):
    base_flux=get_flux_value_tess(teff1)
    b = phoebe.default_binary()
    b.set_value('q',value=q)
    b.set_value('ecc',value=ecc)
    b.set_value('teff', component='primary', value=teff1)
    b.set_value('teff', component='secondary', value=teff2)
    b.set_value('q',value=0.75)#b.add_compute('lc',atm@primary='ck2004',atm@second='ck2004')
#    b.add_dataset('lc',dataset='lc01',ld_func='interp',times = phoebe.linspace(0,1,101),passband='TESS:default',intens_weighting='energy')
    b.add_dataset('lc',dataset='lc01',atm='blackbody',ld_func='interp',times = phoebe.linspace(0,2,101),passband='TESS:default',intens_weighting='energy')
    b.run_compute(irrad_method='none')
#    b['lc01@model'].plot(x='phases',c='blue')
#    b.savefig('lc_kepler.eps')
    times1=b['value@times@lc01@latest@model']
    fluxes1=b['value@fluxes@lc01@latest@model']*base_flux
    return times1,fluxes1

"""
#including foreground bright star
#assumption the distance between binary and us is 1;
#the flux of foreground star 100 times the maximum of flux of binary
#the distance of the binary and foreground star is a function of flux, temperature
"""
def get_foreground_flux_kepler(teff,d):
    base_flux=get_flux_value_kepler(teff)
    d0=1
    #kepler passband flux 
    star_pb1=phoebe.default_star()
    star_pb1.set_value('teff', value=teff)
    star_pb1.add_dataset('lc',dataset='lc01',atm='blackbody',ld_func='interp',times = phoebe.linspace(0,2,101),passband='Kepler:mean',intens_weighting='energy')
#    star_pb1.add_dataset('lc',dataset='lc01',ld_func='interp',times = phoebe.linspace(0,1,101),passband='Kepler:mean',intens_weighting='energy')
    star_pb1.run_compute(irrad_method='none')
    star_pb1_flux=star_pb1['value@fluxes@lc01@latest@model']*(d0/d)**2*base_flux
    return star_pb1_flux

def get_foreground_flux_tess(teff,d):
    base_flux=get_flux_value_tess(teff)
    d0=1
    #kepler passband flux 
    star_pb1=phoebe.default_star()
    star_pb1.set_value('teff', value=teff)
    star_pb1.add_dataset('lc',dataset='lc01',atm='blackbody',ld_func='interp',times = phoebe.linspace(0,2,101),passband='TESS:default',intens_weighting='energy')
#    star_pb1.add_dataset('lc',dataset='lc01',ld_func='interp',times = phoebe.linspace(0,1,101),passband='TESS:default',intens_weighting='energy')
    star_pb1.run_compute(irrad_method='none')
    star_pb1_flux=star_pb1['value@fluxes@lc01@latest@model']*(d0/d)**2*base_flux
    return star_pb1_flux

"""
draw transit light curve
"""
"""
foreground star, solar type teff=6000K, teff1=teff2=6000k,5000K,4000K,7000K
"""

teff1=7000
teff2=7000
teff=6000
d=0.1
d0=1
q=1
ecc=0.0

data1=get_binary_flux_kepler(teff1,teff2,q,ecc)
data2=get_binary_flux_tess(teff1,teff2,q,ecc)
data3=get_foreground_flux_kepler(teff,d)
data4=get_foreground_flux_tess(teff,d)

times=data1[0]
binary_flux_kepler=data1[1]
binary_flux_tess=data2[1]
foreground_flux_kepler=data3
foreground_flux_tess=data4

flux_ratio_k=binary_flux_kepler+foreground_flux_kepler
flux_ratio_k_norm=flux_ratio_k/max(flux_ratio_k)

flux_ratio_t=binary_flux_tess+foreground_flux_tess
flux_ratio_t_norm=flux_ratio_t/max(flux_ratio_t)

fig=plt.figure()
plt.xlabel('Time [d]')
plt.ylabel('flux ratio')
plt.plot(times,flux_ratio_k_norm,c='red',linewidth=0.5,label='kepler')
plt.plot(times,flux_ratio_t_norm,c='blue',linewidth=0.5,label='tess')
plt.legend(loc='lower right')
#plt.yscale('log')
plt.savefig("binary_7k7k_fg_6k.eps")
print(max(abs(flux_ratio_k_norm-flux_ratio_t_norm)))

teff1=6000
teff2=6000
teff=6000
d=0.1
d0=1
q=1
ecc=0.0

data1=get_binary_flux_kepler(teff1,teff2,q,ecc)
data2=get_binary_flux_tess(teff1,teff2,q,ecc)
data3=get_foreground_flux_kepler(teff,d)
data4=get_foreground_flux_tess(teff,d)

times=data1[0]
binary_flux_kepler=data1[1]
binary_flux_tess=data2[1]
foreground_flux_kepler=data3
foreground_flux_tess=data4

flux_ratio_k=binary_flux_kepler+foreground_flux_kepler
flux_ratio_k_norm=flux_ratio_k/max(flux_ratio_k)

flux_ratio_t=binary_flux_tess+foreground_flux_tess
flux_ratio_t_norm=flux_ratio_t/max(flux_ratio_t)

fig=plt.figure()
plt.xlabel('Time [d]')
plt.ylabel('flux ratio')
plt.plot(times,flux_ratio_k_norm,c='red',linewidth=0.5,label='kepler')
plt.plot(times,flux_ratio_t_norm,c='blue',linewidth=0.5,label='tess')
plt.legend(loc='lower right')
#plt.yscale('log')
plt.savefig("binary_6k6k_fg_6k.eps")
print(max(abs(flux_ratio_k_norm-flux_ratio_t_norm)))

teff1=5000
teff2=5000
teff=6000
d=0.1
d0=1
q=1
ecc=0.0

data1=get_binary_flux_kepler(teff1,teff2,q,ecc)
data2=get_binary_flux_tess(teff1,teff2,q,ecc)
data3=get_foreground_flux_kepler(teff,d)
data4=get_foreground_flux_tess(teff,d)

times=data1[0]
binary_flux_kepler=data1[1]
binary_flux_tess=data2[1]
foreground_flux_kepler=data3
foreground_flux_tess=data4

flux_ratio_k=binary_flux_kepler+foreground_flux_kepler
flux_ratio_k_norm=flux_ratio_k/max(flux_ratio_k)

flux_ratio_t=binary_flux_tess+foreground_flux_tess
flux_ratio_t_norm=flux_ratio_t/max(flux_ratio_t)

fig=plt.figure()
plt.xlabel('Time [d]')
plt.ylabel('flux ratio')
plt.plot(times,flux_ratio_k_norm,c='red',linewidth=0.5,label='kepler')
plt.plot(times,flux_ratio_t_norm,c='blue',linewidth=0.5,label='tess')
plt.legend(loc='lower right')
#plt.yscale('log')
plt.savefig("binary_5k5k_fg_6k.eps")
print(max(abs(flux_ratio_k_norm-flux_ratio_t_norm)))


teff1=4000
teff2=4000
teff=6000
d=0.1
d0=1
q=1
ecc=0.0

data1=get_binary_flux_kepler(teff1,teff2,q,ecc)
data2=get_binary_flux_tess(teff1,teff2,q,ecc)
data3=get_foreground_flux_kepler(teff,d)
data4=get_foreground_flux_tess(teff,d)

times=data1[0]
binary_flux_kepler=data1[1]
binary_flux_tess=data2[1]
foreground_flux_kepler=data3
foreground_flux_tess=data4

flux_ratio_k=binary_flux_kepler+foreground_flux_kepler
flux_ratio_k_norm=flux_ratio_k/max(flux_ratio_k)

flux_ratio_t=binary_flux_tess+foreground_flux_tess
flux_ratio_t_norm=flux_ratio_t/max(flux_ratio_t)

fig=plt.figure()
plt.xlabel('Time [d]')
plt.ylabel('flux ratio')
plt.plot(times,flux_ratio_k_norm,c='red',linewidth=0.5,label='kepler')
plt.plot(times,flux_ratio_t_norm,c='blue',linewidth=0.5,label='tess')
plt.legend(loc='lower right')
#plt.yscale('log')
plt.savefig("binary_4k4k_fg_6k.eps")
print(max(abs(flux_ratio_k_norm-flux_ratio_t_norm)))

teff1=6000
teff2=5000
teff=6000
d=0.1
d0=1
q=1
ecc=0.0

data1=get_binary_flux_kepler(teff1,teff2,q,ecc)
data2=get_binary_flux_tess(teff1,teff2,q,ecc)
data3=get_foreground_flux_kepler(teff,d)
data4=get_foreground_flux_tess(teff,d)

times=data1[0]
binary_flux_kepler=data1[1]
binary_flux_tess=data2[1]
foreground_flux_kepler=data3
foreground_flux_tess=data4

flux_ratio_k=binary_flux_kepler+foreground_flux_kepler
flux_ratio_k_norm=flux_ratio_k/max(flux_ratio_k)

flux_ratio_t=binary_flux_tess+foreground_flux_tess
flux_ratio_t_norm=flux_ratio_t/max(flux_ratio_t)

fig=plt.figure()
plt.xlabel('Time [d]')
plt.ylabel('flux ratio')
plt.plot(times,flux_ratio_k_norm,c='red',linewidth=0.5,label='kepler')
plt.plot(times,flux_ratio_t_norm,c='blue',linewidth=0.5,label='tess')
plt.legend(loc='lower right')
#plt.yscale('log')
plt.savefig("binary_6k5k_fg_6k.eps")
print(max(abs(flux_ratio_k_norm-flux_ratio_t_norm)))


teff1=6000
teff2=5000
teff=4000
d=0.1
d0=1
q=1
ecc=0.0

data1=get_binary_flux_kepler(teff1,teff2,q,ecc)
data2=get_binary_flux_tess(teff1,teff2,q,ecc)
data3=get_foreground_flux_kepler(teff,d)
data4=get_foreground_flux_tess(teff,d)

times=data1[0]
binary_flux_kepler=data1[1]
binary_flux_tess=data2[1]
foreground_flux_kepler=data3
foreground_flux_tess=data4

flux_ratio_k=binary_flux_kepler+foreground_flux_kepler
flux_ratio_k_norm=flux_ratio_k/max(flux_ratio_k)

flux_ratio_t=binary_flux_tess+foreground_flux_tess
flux_ratio_t_norm=flux_ratio_t/max(flux_ratio_t)

fig=plt.figure()
plt.xlabel('Time [d]')
plt.ylabel('flux ratio')
plt.plot(times,flux_ratio_k_norm,c='red',linewidth=0.5,label='kepler')
plt.plot(times,flux_ratio_t_norm,c='blue',linewidth=0.5,label='tess')
plt.legend(loc='lower right')
#plt.yscale('log')
plt.savefig("binary_6k5k_fg_4k.eps")
print(max(abs(flux_ratio_k_norm-flux_ratio_t_norm)))


teff1=6000
teff2=5000
teff=7000
d=0.1
d0=1
q=1
ecc=0.0

data1=get_binary_flux_kepler(teff1,teff2,q,ecc)
data2=get_binary_flux_tess(teff1,teff2,q,ecc)
data3=get_foreground_flux_kepler(teff,d)
data4=get_foreground_flux_tess(teff,d)

times=data1[0]
binary_flux_kepler=data1[1]
binary_flux_tess=data2[1]
foreground_flux_kepler=data3
foreground_flux_tess=data4

flux_ratio_k=binary_flux_kepler+foreground_flux_kepler
flux_ratio_k_norm=flux_ratio_k/max(flux_ratio_k)

flux_ratio_t=binary_flux_tess+foreground_flux_tess
flux_ratio_t_norm=flux_ratio_t/max(flux_ratio_t)

fig=plt.figure()
plt.xlabel('Time [d]')
plt.ylabel('flux ratio')
plt.plot(times,flux_ratio_k_norm,c='red',linewidth=0.5,label='kepler')
plt.plot(times,flux_ratio_t_norm,c='blue',linewidth=0.5,label='tess')
plt.legend(loc='lower right')
#plt.yscale('log')
plt.savefig("binary_6k5k_fg_7k.eps")
print(max(abs(flux_ratio_k_norm-flux_ratio_t_norm)))


teff1=7000
teff2=4000
teff=6000
d=0.1
d0=1
q=1
ecc=0.0

data1=get_binary_flux_kepler(teff1,teff2,q,ecc)
data2=get_binary_flux_tess(teff1,teff2,q,ecc)
data3=get_foreground_flux_kepler(teff,d)
data4=get_foreground_flux_tess(teff,d)

times=data1[0]
binary_flux_kepler=data1[1]
binary_flux_tess=data2[1]
foreground_flux_kepler=data3
foreground_flux_tess=data4

flux_ratio_k=binary_flux_kepler+foreground_flux_kepler
flux_ratio_k_norm=flux_ratio_k/max(flux_ratio_k)

flux_ratio_t=binary_flux_tess+foreground_flux_tess
flux_ratio_t_norm=flux_ratio_t/max(flux_ratio_t)

fig=plt.figure()
plt.xlabel('Time [d]')
plt.ylabel('flux ratio')
plt.plot(times,flux_ratio_k_norm,c='red',linewidth=0.5,label='kepler')
plt.plot(times,flux_ratio_t_norm,c='blue',linewidth=0.5,label='tess')
plt.legend(loc='lower right')
#plt.yscale('log')
plt.savefig("binary_7k4k_fg_6k.eps")
print(max(abs(flux_ratio_k_norm-flux_ratio_t_norm)))


teff1=7000
teff2=4000
teff=20000
d=0.1
d0=1
q=1
ecc=0.0

data1=get_binary_flux_kepler(teff1,teff2,q,ecc)
data2=get_binary_flux_tess(teff1,teff2,q,ecc)
data3=get_foreground_flux_kepler(teff,d)
data4=get_foreground_flux_tess(teff,d)

times=data1[0]
binary_flux_kepler=data1[1]
binary_flux_tess=data2[1]
foreground_flux_kepler=data3
foreground_flux_tess=data4

flux_ratio_k=binary_flux_kepler+foreground_flux_kepler
flux_ratio_k_norm=flux_ratio_k/max(flux_ratio_k)

flux_ratio_t=binary_flux_tess+foreground_flux_tess
flux_ratio_t_norm=flux_ratio_t/max(flux_ratio_t)

fig=plt.figure()
plt.xlabel('Time [d]')
plt.ylabel('flux ratio')
plt.plot(times,flux_ratio_k_norm,c='red',linewidth=0.5,label='kepler')
plt.plot(times,flux_ratio_t_norm,c='blue',linewidth=0.5,label='tess')
plt.legend(loc='lower right')
#plt.yscale('log')
plt.savefig("binary_7k4k_fg_20k.eps")
print(max(abs(flux_ratio_k_norm-flux_ratio_t_norm)))

teff1=7000
teff2=4000
teff=10000
d=0.1
d0=1
q=1
ecc=0.0

data1=get_binary_flux_kepler(teff1,teff2,q,ecc)
data2=get_binary_flux_tess(teff1,teff2,q,ecc)
data3=get_foreground_flux_kepler(teff,d)
data4=get_foreground_flux_tess(teff,d)

times=data1[0]
binary_flux_kepler=data1[1]
binary_flux_tess=data2[1]
foreground_flux_kepler=data3
foreground_flux_tess=data4

flux_ratio_k=binary_flux_kepler+foreground_flux_kepler
flux_ratio_k_norm=flux_ratio_k/max(flux_ratio_k)

flux_ratio_t=binary_flux_tess+foreground_flux_tess
flux_ratio_t_norm=flux_ratio_t/max(flux_ratio_t)

fig=plt.figure()
plt.xlabel('Time [d]')
plt.ylabel('flux ratio')
plt.plot(times,flux_ratio_k_norm,c='red',linewidth=0.5,label='kepler')
plt.plot(times,flux_ratio_t_norm,c='blue',linewidth=0.5,label='tess')
plt.legend(loc='lower right')
#plt.yscale('log')
plt.savefig("binary_7k4k_fg_10k.eps")
print(max(abs(flux_ratio_k_norm-flux_ratio_t_norm)))

teff1=10000
teff2=4000
teff=6000
d=0.1
d0=1
q=1
ecc=0.0

data1=get_binary_flux_kepler(teff1,teff2,q,ecc)
data2=get_binary_flux_tess(teff1,teff2,q,ecc)
data3=get_foreground_flux_kepler(teff,d)
data4=get_foreground_flux_tess(teff,d)

times=data1[0]
binary_flux_kepler=data1[1]
binary_flux_tess=data2[1]
foreground_flux_kepler=data3
foreground_flux_tess=data4

flux_ratio_k=binary_flux_kepler+foreground_flux_kepler
flux_ratio_k_norm=flux_ratio_k/max(flux_ratio_k)

flux_ratio_t=binary_flux_tess+foreground_flux_tess
flux_ratio_t_norm=flux_ratio_t/max(flux_ratio_t)

fig=plt.figure()
plt.xlabel('Time [d]')
plt.ylabel('flux ratio')
plt.plot(times,flux_ratio_k_norm,c='red',linewidth=0.5,label='kepler')
plt.plot(times,flux_ratio_t_norm,c='blue',linewidth=0.5,label='tess')
plt.legend(loc='lower right')
#plt.yscale('log')
plt.savefig("binary_10k4k_fg_6k.eps")
print(max(abs(flux_ratio_k_norm-flux_ratio_t_norm)))
