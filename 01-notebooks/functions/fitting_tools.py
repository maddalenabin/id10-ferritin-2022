from scipy.optimize import curve_fit
import numpy as np

def fit(function,x,y,p0=None,sigma=None,bounds=[None,None]):
    """Fits a function and return the fit resulting parameters and curve
    """
    popt,pcov = curve_fit(function,x,y,p0=p0,sigma=sigma,bounds=bounds)
    xc = np.linspace(min(x),max(x),10000)
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


def exponential(x, beta, tau):
    """ Single exponential (brownian motion) """
    return np.abs(beta) * np.exp( -2*x/(np.abs(tau)) )

def exponential_kww(x, beta, tau, kww):
    """ Stretched/compressed exponential """
    return np.abs(beta) * np.exp( -2*(x/(np.abs(tau)))**kww )

def linear(x, m):
    """ Linear function """
    return m*x


def gaussian(x, a, mean, sigma):
    """ Gaussian """
    return a * np.exp( ((x-mean)/sigma)**2 )

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

