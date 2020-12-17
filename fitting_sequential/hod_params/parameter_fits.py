import numpy as np
import HOD_curves as hod
from scipy.optimize import curve_fit


def fit_M(magnitude, logM, mask=None, function=hod.M_function, p0=None):
    """
    Fit smooth curve to Mmin or M1
    
    magnitude: array of absolute magnitudes
    logM:      array of logMmin or logM1
    mask:      optional boolean array, for removing samples from the fit
    function:  function to fit
    p0:        initial parameter guess
    """
    if mask is None:
        popt, pcov = curve_fit(function, logM, magnitude, p0=p0)
    else:
        popt, pcov = curve_fit(function, logM[mask], magnitude[mask], p0=p0)
    print("popt", popt)
    
    # get new values of Mmin or M1
    logM_new = np.zeros(len(magnitude))
    logM_fine = np.arange(10,17,0.00001)
    mags_fine = function(logM_fine, *popt)
    
    for i in range(len(magnitude)):
        logM_new[i] = logM_fine[mags_fine <= magnitude[i]][0]
        
    return logM_new


def fit_sigma(magnitude, sigma, mask=None, function=hod.sigma_function, p0=None):
    """
    Fit smooth curve to sigma
    
    magnitude: array of absolute magnitudes
    logM:      array of sigma
    mask:      optional boolean array, for removing samples from the fit
    function:  function to fit
    p0:        initial parameter guess
    """
    
    if mask is None:
        popt, pcov = curve_fit(function, magnitude, sigma, p0=p0)
    else:
        popt, pcov = curve_fit(function, magnitude[mask], sigma[mask], p0=p0)
    print("popt", popt)
    
    sigma_new = function(magnitude, *popt)
    
    return sigma_new


def fit_M0(magnitude, logM0, mask=None, function=hod.M0_function, p0=None):
    """
    Fit smooth curve to M0
    
    magnitude: array of absolute magnitudes
    logM:      array of logM0
    mask:      optional boolean array, for removing samples from the fit
    function:  function to fit
    p0:        initial parameter guess
    """
    
    log_function = lambda x,a,b,c,d: np.log10(function(x,a,b,c,d))
    
    if mask is None:
        popt, pcov = curve_fit(log_function, magnitude, logM0, p0=p0)
    else:
        popt, pcov = curve_fit(log_function, magnitude[mask], logM0[mask], p0=p0)
    print("popt", popt)
        
    logM0_new = log_function(magnitude, *popt)
    
    return logM0_new


def fit_alpha(magnitude, alpha, mask=None, function=hod.alpha_function, p0=None):
    """
    Fit smooth curve to alpha
    
    magnitude: array of absolute magnitudes
    logM:      array of alpha
    mask:      optional boolean array, for removing samples from the fit
    function:  function to fit
    p0:        initial parameter guess
    """
    
    if mask is None:
        popt, pcov = curve_fit(function, magnitude, alpha, p0=p0)
    else:
        popt, pcov = curve_fit(function, magnitude[mask], alpha[mask], p0=p0)
    print("popt", popt)
    
    alpha_new = function(magnitude, *popt)
    
    return alpha_new
