import numpy as np
from scipy.optimize import curve_fit

def fit(function,x,y,p0=None,sigma=None,bounds=[None,None]):
    """Fits a function and return the fit resulting parameters and curve
    """
    popt,pcov = curve_fit(function,x,y,p0=p0,sigma=sigma,bounds=bounds)
    xc = np.linspace(min(x),max(x),10000)
    curve = function(xc,*popt)
    perr = np.sqrt(np.diag(pcov))
    return popt,xc,curve,perr