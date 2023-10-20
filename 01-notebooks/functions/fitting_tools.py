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


def exponential(x, beta, tau):
    """ Single exponential (brownian motion) """

    return np.abs(beta) * np.exp( -2*x/(np.abs(tau)) )

def exponential_kww(x, beta, tau, kww):
    """ Stretched/compressed exponential """
    
    return np.abs(beta) * np.exp( -2*(x/(np.abs(tau)))**kww )

def linear(x, m):
    """ Linear function """
    
    return m*x


def gaussian(x, a, mean, sigma, tau):
    """ Gaussian """

    return a * np.exp( ((x-mean)/sigma)**2 )


# add Arrhenius