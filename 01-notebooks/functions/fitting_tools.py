import numpy as np
import pandas as pd
import pickle
from scipy.optimize import curve_fit

def fit(function,x,y,p0=None,sigma=None,bounds=[None,None]):
    """Fits a function and return the fit resulting parameters and curve
    """
    popt,pcov = curve_fit(function,x,y,p0=p0,sigma=sigma,bounds=bounds)
    xc = np.linspace(min(x),max(x),10000)
    curve = function(xc,*popt)
    perr = np.sqrt(np.diag(pcov))
    return popt,xc,curve,perr


def D_coeff(T, eta, Rh=6.6e-9):
    """Calculate the diffusion coefficient in nm^2/us from Stokes-Einsten equation
    Args
        T: temperature in K
        eta: viscosity of the solvent in N*s/m^2 
        Rh: hydrodynamic radius in m. Default value is for ferritin
    """
    if T == 243:
        T = 243
    kb = 1.380649e-23 # m2 kg s-2 K-1
    D = kb*T / (6*np.pi*eta*Rh) # m^2 / s
    
    return D*1e12


def eta_50v_gly():
    """Viscosity values in Ns/m2 of water/glycerol mixture with 50v% glycerol.
       http://www.met.reading.ac.uk/~sws04cdw/viscosity_calc.html
    """
    viscosity_50v_gly = {270: 0.026578, 260: 0.052163, 250: 0.11934, 240: 0.33719, 230: 1.2893, 220: 7.7852, 210: 98.826}

    return viscosity_50v_gly

def eta_60v_gly():
    """Viscosity values in Ns/m2 of water/glycerol mixture with 60v% glycerol.
       http://www.met.reading.ac.uk/~sws04cdw/viscosity_calc.html
    """
    viscosity_60v_gly = {270: 0.061359, 260: 0.13400, 250: 0.34857, 240: 1.1499, 230: 5.3002, 220: 40.015, 210: 657.37}

    return viscosity_60v_gly


def update_D_coeff(T, Tr, D, dD, c, run):
    """Updates the diffusion coefficient calculated from the data
    """
    DiffC = pd.read_pickle("/cfs/home/mabi3848/id10-ferritin-2022/03-source/diffusion_coefficient")
    exists = DiffC['run'].str.contains(run)
    
    if exists:
        pass
    else:
        new_row = {'temperature': T, 'transmission': Tr, 'D': D, 'dD': dD, 'c': c, 'run': run}
        DiffC2 = DiffC.append(new_row, ignore_index=True)
        DiffC2.to_pickle("/cfs/home/mabi3848/id10-ferritin-2022/03-source/diffusion_coefficient")

    return