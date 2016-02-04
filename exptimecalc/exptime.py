"""
This is an exposure time calculator for Apache Point Observatory (APO).
Copyright 2016 S. Drew Chojnowski (NMSU).
"""
import math
import numpy as np

# Vega flux zeropoints
# borrowed from http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
# vega_*: (erg/s/cm^2/angstrom)
vega_u=4.175e-9     
vega_b=6.32e-9      
vega_v=3.631e-9     
vega_r=2.177e-9     
vega_i=1.126e-9     
vega_j=0.3147e-9    
vega_h=0.1138e-9    
vega_k=0.03961e-9

# Johnson filters
# borrowed from http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
# units=angstroms
wave_u=3600.0  ; width_u=600.0
wave_b=4380.0  ; width_b=900.0
wave_v=5450.0  ; width_v=850.0
wave_r=6410.0  ; width_r=1500.0
wave_i=7980.0  ; width_i=1500.0
wave_j=12200.0  ; width_j=2600.0
wave_h=16300.0  ; width_h=2900.0
wave_k=21900.0  ; width_k=4100.0

# constants
hplanck=6.6260755e-27 # erg s
cspeed=2.99792458e18  # angstrom/s

def exptime(telescope='nmsu1m',instrument=None,exptime=1,wavelength=5500,bandwidth=1000,starmag=22,band='V',
            magsystem='Vega',seeing=1.0,airmass=1.0,lunarphase=7.0,silent=False):

    if magsystem is None:
        choice=raw_input('Enter 1 for Vega mag. system, 2 for AB mag. system: ')
        if choice=='1': magsystem='Vega'
        if choice=='2': magsystem='AB'

    sys_efficiency=get_system_efficiency('nmsu1m')
    atm_trans=get_atmospheric_transmission(wavelength=wavelength)
    if telescope=='nmsu1m': mirror_area=get_mirror_area(1,'m')
    if telescope=='3.5m': mirror_area=get_mirror_area(3.5,'m')
    flux=get_flux(starmag=starmag,band=band,magsystem=magsystem)
    int_lambda=0.5*((wavelength+bandwidth/2.)**2-(wavelength-bandwidth/2.)**2)
    signal=mirror_area*exptime*atm_trans*sys_efficiency*flux*(1./hplanck)*(1./cspeed)*int_lambda

    if silent is False:
        print '===================================================================='  
        print '# of photons collected in a ',exptime,' second observation'
        print 'of a ',band,'=',starmag,' object by the ',telescope,' telescope:\n'
        print signal
        print '\nSNR = ',signal/math.sqrt(signal)

    return signal



def get_system_efficiency(telescope=None,instrument=None):
    if telescope=='nmsu1m': se=0.5
    if telescope=='apo3.5m': se=0.5

    return se

def get_atmospheric_transmission(wavelength=None):
    at=0.8

    return at

def get_mirror_area(msize=None,munits=None):
    if munits=='m': msize=msize*100.0
    if munits=='cm': msize=msize
    ma=math.pi*((msize/2.)**2)

    return ma

def get_flux(wavelength=None,starmag=None,band=None,magsystem=None):
    if magsystem=='Vega':
        if band=='U': flux=vega_u*(10**(starmag/(-2.5)))
        if band=='B': flux=vega_b*(10**(starmag/(-2.5)))
        if band=='V': flux=vega_v*(10**(starmag/(-2.5)))
        if band=='R': flux=vega_r*(10**(starmag/(-2.5)))
        if band=='I': flux=vega_i*(10**(starmag/(-2.5)))
        if band=='J': flux=vega_j*(10**(starmag/(-2.5)))
        if band=='H': flux=vega_h*(10**(starmag/(-2.5)))
        if band=='K': flux=vega_k*(10**(starmag/(-2.5)))

    return flux







