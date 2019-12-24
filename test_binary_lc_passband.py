#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 16:18:03 2019

@author: dyz
"""

import phoebe
from phoebe import u
import numpy as np
import matplotlib.pyplot as plt

logger = phoebe.logger(clevel='WARNING')

kepler = phoebe.get_passband('Kepler:mean')
tess = phoebe.get_passband('TESS:default')

#calculate the two passband flux
teffs = np.linspace(3500, 8000, 100)
fig=plt.figure()
plt.xlabel('Temperature [K]')
plt.ylabel('Inorm [W/m^2/A]')
plt.plot(teffs, kepler.Inorm(teffs, atm='blackbody', ld_func='linear', ld_coeffs=[0.0]), label='Kepler')
plt.plot(teffs, tess.Inorm(teffs, atm='blackbody', ld_func='linear', ld_coeffs=[0.0]), label='TESS')
plt.legend(loc='lower right')

plt.savefig("two_passband_flux.eps")
#plt.show()

"""
b = phoebe.default_binary()
times = np.linspace(0,1,21)
b.add_dataset('lc', times=times, dataset='lc01')
b.add_dataset('rv', times=times, dataset='rv01')
b.add_dataset('mesh', times=times, columns=['visibilities', 'intensities@lc01', 'rvs@rv01'], dataset='mesh01')
b.run_compute(irrad_method='none')
b['lc01@model'].plot(axpos=221)
b['rv01@model'].plot(c={'primary': 'blue', 'secondary': 'red'}, linestyle='solid', axpos=222)
b['mesh@model'].plot(fc='intensities@lc01', ec='None', axpos=425)
b['mesh@model'].plot(fc='rvs@rv01', ec='None', axpos=427)
b['mesh@model'].plot(fc='visibilities', ec='None', y='ws', axpos=224)

fig = plt.figure(figsize=(11,4))
b.savefig('animation_binary_complete.gif', fig=fig, tight_layouot=True, draw_sidebars=False, animate=True, save_kwargs={'writer': 'imagemagick'})
plt.show()
"""


b = phoebe.default_binary()
#b['value@teff@primary@component']=6000
#b['value@teff@secondary@component']=5000

b.set_value('teff', component='secondary', value=5000)
b.set_value('q',value=0.75)#b.add_compute('lc',atm@primary='ck2004',atm@second='ck2004')
b.add_dataset('lc',dataset='lc01',ld_func='interp',times = phoebe.linspace(0,2,101),passband='Kepler:mean',intens_weighting='energy')
b.run_compute(irrad_method='none')
b['lc01@model'].plot(x='phases',c='blue')
b.savefig('lc_kepler.eps')
times1=b['value@times@lc01@latest@model']
fluxes1=b['value@fluxes@lc01@latest@model']
#plt.legend(loc='lower right')

c = phoebe.default_binary()
c.set_value('teff', component='secondary', value=5000)
c.set_value('q',value=0.75)#b.add_compute('lc',atm@primary='ck2004',atm@second='ck2004')
c.add_dataset('lc',dataset='lc01',ld_func='interp',times = phoebe.linspace(0,2,101),passband='TESS:default',intens_weighting='energy')
c.run_compute(irrad_method='none')
c['lc01@model'].plot(x='phases',c='red')
c.savefig('lc_TESS.eps')
times2=c['value@times@lc01@latest@model']
fluxes2=c['value@fluxes@lc01@latest@model']


fig=plt.figure()
plt.xlabel('Time [d]')
plt.ylabel('flux')
plt.plot(times1,fluxes1,c='red',linewidth=0.5,label='kepler')
plt.plot(times2,fluxes2,c='blue',linewidth=0.5,label='tess')
plt.legend(loc='lower right')
plt.savefig("kepler_tess_cl.eps")

#including foreground bright star
#assumption the distance between binary and us is 1;
#the flux of foreground star 100 times the maximum of flux of binary
#the distance of the binary and foreground star is a function of flux, temperature

foreground_t=6000
d0=1
d=0.1
#kepler passband flux 
star_pb1=phoebe.default_star()
star_pb1.set_value('teff', value=6000)
star_pb1.add_dataset('lc',dataset='lc01',ld_func='interp',times = phoebe.linspace(0,2,101),passband='Kepler:mean',intens_weighting='energy')
star_pb1.run_compute(irrad_method='none')
star_pb1_flux=star_pb1['value@fluxes@lc01@latest@model']*(d0/d)**2

#Tess passband flux
star_pb2=phoebe.default_star()
star_pb2.set_value('teff', value=6000)
star_pb2.add_dataset('lc',dataset='lc01',ld_func='interp',times = phoebe.linspace(0,2,101),passband='TESS:default',intens_weighting='energy')
star_pb2.run_compute(irrad_method='none')
star_pb2_flux=star_pb1['value@fluxes@lc01@latest@model']*(d0/d)**2



fig=plt.figure()
plt.xlabel('Time [d]')
plt.ylabel('flux ratio')
plt.plot(times1,fluxes1/star_pb1_flux,c='red',linewidth=0.5,label='kepler')
plt.plot(times2,fluxes2/star_pb2_flux,c='blue',linewidth=0.5,label='tess')
plt.legend(loc='lower right')
plt.yscale('log')
plt.savefig("kepler_tess_cl_2.eps")


#rv lc
a=phoebe.default_binary()
a.set_value('teff', component='secondary', value=5000)
a.set_value('q',value=0.75)
a.add_dataset('lc', passband='Kepler:mean',times = phoebe.linspace(0,2,101))
a.add_dataset('rv', passband='Kepler:mean',times= phoebe.linspace(0,2,101))
a.run_compute()
a.plot(x='phases')
a.savefig('cl_rv_kepler.eps')


"""
times = np.linspace(0,2,101)
a.add_dataset('mesh', times=times, columns=['visibilities', 'intensities@lc01', 'rvs@rv01'], dataset='mesh01')
a.run_compute()
a['lc01@model'].plot(axpos=221)
a['rv01@model'].plot(c={'primary': 'blue', 'secondary': 'red'}, linestyle='solid', axpos=222)
a['mesh@model'].plot(fc='intensities@lc01', ec='None', axpos=425)
a['mesh@model'].plot(fc='rvs@rv01', ec='None', axpos=427)
a['mesh@model'].plot(fc='visibilities', ec='None', y='ws', axpos=224)
fig = plt.figure(figsize=(11,4))
a.savefig('cl_rv_kepler.eps', fig=fig, tight_layouot=True, draw_sidebars=False, animate=True, save_kwargs={'writer': 'imagemagick'})
"""




"""
a = phoebe.default_binary()
a['value@teff@primary@component']=6000
a['value@teff@secondary@component']=5000

a.add_dataset('lc',dataset='lc01',ld_func='interp',passband='Kepler:mean',intens_weighting='energy')#,l3=10000)
#a.add_compute('lc',model='Kepler')
a.set_value('times',np.linspace(0,1,21))
a.run_compute(irrad_method='none')

a.add_dataset('lc',dataset='lc02',ld_func='interp',passband='TESS:default',intens_weighting='energy')#,l3=10000)
#a.set_value('times',np.linspace(0,1,21))
#a.add_compute('lc',model='TESS')

a.run_compute(irrad_method='none')
a['lc01@model'].plot(c='blue')
a['lc02@model'].plot(c='red')
a.savefig('lc_two_passband.eps')
"""

"""
def kepler_passband_cl(t_primary,t_secondary,ecc,background_flux):
    b = phoebe.default_binary()
#    b['q'] = q
    b['ecc'] = ecc
    b['value@teff@primary@component']=t_primary
    b['value@teff@secondary@component']=t_primary    
    b.add_dataset('lc',dataset='lc01',ld_func='interp',passband='Kepler:mean',intens_weighting='energy')
    #b.add_compute('lc',atm@primary='ck2004',atm@second='ck2004')
    b.set_value('times',np.linspace(0,1,101))
    b.run_compute(irrad_method='none')
    times=b['value@times@lc01@latest@model']
    fluxes=b['value@fluxes@lc01@latest@model']
    return times,fluxes

def tess_passband_cl(t_primary,t_secondary,ecc,background_flux):
    b = phoebe.default_binary()
#    b['q'] = q
    b['ecc'] = ecc
    b['value@teff@primary@component']=t_primary
    b['value@teff@secondary@component']=t_primary    
    b.add_dataset('lc',dataset='lc01',ld_func='interp',passband='TESS:default',intens_weighting='energy')
    #b.add_compute('lc',atm@primary='ck2004',atm@second='ck2004')
    b.set_value('times',np.linspace(0,1,101))
    b.run_compute(irrad_method='none')
    times=b['value@times@lc01@latest@model']
    fluxes=b['value@fluxes@lc01@latest@model']
    return times,fluxes


#q=1.0
ecc=0.0
t_primary=7000
t_secondary=4000
back_flux=0


kepler_times,kepler_fluxes=kepler_passband_cl(t_primary,t_secondary,ecc,back_flux)
kepler_flux_ratio=(kepler_fluxes-back_flux)/kepler_fluxes

tess_times,tess_fluxes=tess_passband_cl(t_primary,t_secondary,ecc,back_flux)
tess_flux_ratio=(tess_fluxes-back_flux)/tess_fluxes

fig=plt.figure()
plt.xlabel('Time [d]')
plt.ylabel('flux')
plt.plot(kepler_times,kepler_fluxes,c='blue',linewidth=0.5,label='kepler band')
plt.plot(tess_times,tess_fluxes,c='red',linewidth=0.5,label='tess band')
plt.legend(loc='lower right')
plt.savefig("two_passband_cl.eps")
"""