##################################################################
# lightunits(wavelength,frequency,energy)
#
# This function converts between wavelength, frequency, and energy.
# If a wavelength is provided, frequency and energy are calculated,
# and so on.
#
# Input/output units:
#     WAVELENGTH: angstroms (10^-10 m)
#     FREQUENCY: GHz (10^9 Hz)
#     ENERGY: eV 
##################################################################

import math
import numpy as np
import sys
import scipy

def lightunits(wavelength=0,frequency=0,energy=0):
    cspeed=2.99792458e8 # m/s
    hplanck=4.135667662e-15 # eV s

    # find out which value has been provided
    vals=np.array([wavelength,frequency,energy])
    gdvals=np.where(vals!=0)

    # convert values
    if gdvals[0][0]==0:
        print 'frequency: ',(cspeed/(wavelength*(1.0e-10)))/1.0e9,' (GHz)'
        print 'energy: ',(hplanck*cspeed)/(wavelength*(1.0e-10)),' (eV)'
    if gdvals[0][0]==1:
        print 'wavelength: ',(cspeed/(frequency*(1.0e9)))/1.0e-10,' (Angstroms)'
        print 'energy: ',hplanck*(frequency*(1.0e9)),' (eV)'
    if gdvals[0][0]==2:
        print 'wavelength: ',((cspeed*hplanck)/energy)/1.0e-10,' (Angstroms)'
        print 'frequency: ',(energy/hplanck)/1.0e9,' (GHz)'

    return


##################################################################
# fluxunits(wavelength,flambda,fnu)
#
# This function converts between flux per unit wavelength to flux
# per unit frequency, and vice versa, for a given wavelength.
#
# Input/output units:
#     WAVELENGTH: angstroms (10^-10 m)
#     FLAMBDA: erg/cm^2/s/angstrom
#     FNU: erg/cm^2/s/Hz
##################################################################

def fluxunits(wavelength=3000,f_lambda=0.0,f_nu=0):
    cspeed=2.99792458e18 # angstroms/s

    if f_lambda==0 and f_nu==0:
        print 'No flux value provided, so defaulting to f_lambda=1 @ 3000 angstroms'
        f_lambda=1.0

    if f_lambda==0:
        print 'F_lambda:',f_nu*(cspeed/(wavelength**2)),' erg/cm^2/s/angstrom'

    if f_nu==0:
        print 'F_nu:',f_lambda*((wavelength**2)/cspeed),' erg/cm^2/s/Hz'

    return

##################################################################
# abnu2stmag(m_lambda,m_nu)
#
# This function converts between monochromatic apparent magnitudes 
# in the # ABNU (frequency) and STMAG (wavelength) systems.
#
# Input/output units:
#     WAVELENGTH: angstroms (10^-10 m)
#     M_LAMBDA: mags
#     M_NU: mags
##################################################################

def abnu2stmag(wavelength=5500,m_lambda=0,m_nu=0):
    if m_lambda==0 and m_nu==0:
        print 'No mag value provided, so defaulting to m_lambda=1 @ 5500 angstroms'
        f_lambda=1.0

    if m_nu==0:
        print 'm_nu: ',m_lambda-(5.0*math.log10(wavelength))+18.70175,' mags'

    if m_lambda==0:
        print 'm_lambda: ',m_nu+(5.0*math.log10(wavelength))-18.70175,' mags'

    return




