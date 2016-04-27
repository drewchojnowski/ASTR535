import numpy as np
import math
from astropy.io import ascii
from astropy import constants as const
from astropy import units as u
import matplotlib.pyplot as plt

vegaspectrum='../vega_c95.txt'

def test():
#    instr=get_instr()
#    mag,airmass,teff,moonphase,filt,w1,w2,wcent=get_objinfo()

    wcent=6250.0
    teff=5000.0
    mag=9.0
    w1=5500.0; w2=7000.0
    vega = ascii.read(vegaspectrum)
    w=np.array(vega['wavelength'])
    f=np.array(vega['flux'])
    gd=np.where((w>wcent-0.6) & (w<wcent+0.6))

    refflux=f[gd][0]
    wave,objflux=planck(teff=teff,refflux=refflux,mag=mag,w1=w1,w2=w2,wcent=wcent)

    return wave,objflux

def get_instr():
    print('Which instrument are you using?')
    print('(1) DIS')
    print('(2) ARCTIC')
    a=raw_input('')
    if a=='1': instr='dis'
    if a=='2': instr='arctic'

    return instr

def get_objinfo():
    # note: need to error check the user input
    print('What is the magnitude of the object?')
    mag=raw_input('')
    print('What is the airmass?')
    airmass=raw_input('')
    print('What is the expected T_eff?')
    teff=raw_input('')
    print('What is the lunar phase?')
    moonphase=raw_input('')
    print('What filter are you using?')
    filt=raw_input('')
    print('What grating are you using?')
    grating=raw_input('')

    mag=9.0
    airmass=1.0
    teff=5000.0
    moonphase=0.0
    filt='sdss_r'
    w1=5500.0
    w2=7000.0
    wcent=mean([w1,w2])
    
    return mag,airmass,teff,moonphase,filt,w1,w2,wcent

'''
       The reference spectrum will be Vega, as taken from a table of values 
       and formatted in an array for this program.
'''
h = const.h.cgs.value
c = const.c.cgs.value
k = const.k_B.cgs.value
bandpasses = {'U':(3000,4500,3656), 'B':(3500,6000,4353), 'V':(4500,7500,5477), 
           'R':(5000,9500,6349),'I':(6500,12500,8797)}

def planck(teff,refflux,mag,band=None,w1=None,w2=None,wcent=None):
    '''
       Inputs:  Target's effective temperature and magnitude at a specific wavelength, the 
                desired bandpass, and a reference flux.
       Ouputs:  Array of wavelengths within the bandpass and the corresponding fluxes.
    '''
    abscenter = refflux*10**(-0.4*mag)
    if band is None: 
        lolim,hilim,center=w1,w2,wcent
    else:
        lolim,hilim,center = bandpasses[band]
    wave = np.arange(lolim,hilim)
    modwave = wave*1e-8
    unflux = (2.0*h*c**2/modwave**5) * (np.exp(h*c/(modwave*k*teff))-1.0)**(-1)

    pairs = dict(zip(wave,unflux))                # combines two lists into dictionary
    uncenter = pairs[center]                      # unnormalized flux at band's central wavelength
    N = abscenter/uncenter                        # normalization constant
    objflux = N*unflux                           # normalized flux curve
    return wave,objflux

def planck_graph():
    x,y = planck(5000,3.6e-9,0,'U')
    plt.plot(x,y)
    plt.show()


# Module to calculate q (system efficiency)

# Open files from which to read in system data

def mirrorThru(waveR,reflect):
    '''
    Calculates reflectivity of mirrors as a function of wavelength.
    Returns: Transmission of mirrors
    '''
    r = open('../Input_files/al_reflectivity.txt')
    waveR,reflect = np.loadtxt(r, unpack=True)
    mirrors = (interpolate(waveR,reflect))**3.0
    # To the third power - for each of the three mirrors on the 3.5m
    return mirrors

def filThru(waveF,filt):
    '''
    Returns: filter transmission as a function of wavelength.
    '''
    f = open('filter.txt')
    waveF,filt = np.loadtxt(f,unpack=True)
    filt = interpolate(waveF,filt)
    return filt

def detecThru(waveQE,qe):
    '''
    Calculates detector transmission as a function of wavelength.
    Takes into consideration: Quantum efficiency, readout noise,
    electron gain, plate scale, seeing.
    '''
    plate = 0.228  # Plate scale for ARCTIC (arcsec/pix) (2x2 binning)
    gain = 2.00    # electrons/DN
    rn = 3.7       # electrons
    qeff = interpolate(waveQE,qe)
    detec = ???
    return detec

def calc_q(mirrors,filt,detec):
    '''
    Multiplies all factors of system efficiency.
    Returns: q as a function of wavelength.
    '''
    qe = open('qeff.txt')
    waveQE,qe = np.loadtxt(f,unpack=True)
    net_q = mirrors*filt*detec
    return net_q






