from scipy.optimize import curve_fit
import numpy as np

def fit(function,x,y,p0=None,sigma=None,bounds=[None,None]):
    """Fits a function and return the fit resulting parameters and curve
    """
    popt,pcov = curve_fit(function,x,y,p0=p0,sigma=sigma,bounds=bounds)
    xc = np.linspace(min(x),max(x),1000)
    curve = function(xc,*popt)
    perr = np.sqrt(np.diag(pcov))
    return popt,xc,curve,perr

def fit2(function,x,y,p0=None,sigma=None,bounds=[None,None]):
    """Fits a function and return the fit resulting parameters and curve
    """
    popt,pcov = curve_fit(function,x,y,p0=p0,sigma=sigma,bounds=bounds)
    xc = np.linspace(min(x),max(x),len(x))
    curve = function(xc,*popt)
    perr = np.sqrt(np.diag(pcov))
    return popt,xc,curve,perr

def fit3(function,x,y,p0=None,sigma=None,bounds=[None,None]):
    """Fits a function and return the fit resulting parameters and curve
    """
    popt,pcov = curve_fit(function,x,y,p0=p0,sigma=sigma,bounds=bounds)
    # xc = np.linspace(min(x),max(x),len(x))
    curve = function(x,*popt)
    perr = np.sqrt(np.diag(pcov))
    return popt,x,curve,perr


def exponential(x, beta, tau):
    """ Single exponential (brownian motion) """
    return np.abs(beta) * np.exp( -2*x/(np.abs(tau)) )

def exponential_kww(x, beta, tau, kww):
    """ Stretched/compressed exponential """
    return np.abs(beta) * np.exp( -2*(x/(np.abs(tau)))**kww )

# def exponential_kww(x, beta, tau, kww):
#     """ Stretched/compressed exponential with fixed contrast"""
#     return np.abs(beta) * np.exp( -2*(x/(np.abs(tau)))**kww )

def linear(x, m):
    """ Linear function """
    return m*x

def linear_q(x,m,q):
    """ Linear function with offset """
    return m * x + q

def gaussian(x, a, mean, sigma):
    """ Gaussian """
    return a * np.exp( ((x-mean)/sigma)**2 )

def gaussian_off(x, a, mean, sigma, c):
    """ Gaussian with baseline offset"""
    return a * np.exp( -((x-mean)/sigma)**2 ) + c

def arrhenius(x,t1,Ea):
    """  Arrhenius function """
    Kb = 1#1.380649e-23 #m2 kg s-2 K-1 = J/K
    # Ea = popt2[1] * Kb / 1e3 * Na # if you use Kb = 1, then here use the real Kb, where Na = 6.02214076e23
    return t1 * np.exp( Ea/(x*Kb) )

def VFT(x, a, D, T0):
    """ Vogel–Fulcher–Tammann for fragile glass formers 
    Args: a amplitude, D fragility parameter, T0 ideal glass transition temperature
    """
    return a * np.exp(D * T0 / (x-T0))

def VFT2(x, a, B, T0):
    """ Vogel–Fulcher–Tammann for fragile glass formers 
    Args: a amplitude, D fragility parameter, T0 ideal glass transition temperature
    https://link.aps.org/doi/10.1103/PhysRevLett.96.145502
    """
    return a * np.sqrt(x) * np.exp( B / (x-T0))

def eyring(T, H, S):
    """ Eyring equation
    Args: T temperature, H enthalpy, S entropy
    """
    kb = 1.380649e-23 # m2 kg s-2 K-1
    R = 8.314 # J / mol·K # gas constant
    h =  6.626e-34 # Plank's constant  kg⋅m2⋅s−1 = J⋅s
    return H / R / T + S/R + np.log(kb/h)

def viscosity_predic_fit(T, eta0, alpha, Tg):
    " this is valid for w=0.5 from Gonzalez et al J. Chem. Eng. Data 2011, 56, 1397–1406"
    eta = eta0 * np.exp( (28.75 - np.log(eta0)) * (Tg / T)**alpha ) # 
    return eta