import numpy as np
from scipy.interpolate import splrep, splev


def logMmin_to_sigma(magnitude, logMmin, filename):
    """
    Use degeneracy to convert logMmin to sigma_logM
    """
    sigma_new = np.zeros(len(magnitude))
    
    arr = np.loadtxt(filename)
    
    for i in range(len(magnitude)):
        sig = arr[:,0]
        logM = arr[:,i+1]
        try:
            tck = splrep(logM, sig)
            sigma_new[i] = splev(logMmin[i], tck)
        except:
            # for the faint samples, logMmin can be larger than the maximum value along the degeneracy curve
            # in this case, just find where the degeneracy curve reaches a maximum
            diff = np.absolute(logMmin[i] - logM)
            idx = np.argmin(diff)
            sigma_new[i] = sig[idx]
            
    return sigma_new



def sigma_to_logMmin(magnitude, sigma, filename):
    """
    Use degeneracy to convert sigma_logM to logMmin
    """
    logMmin_new = np.zeros(len(magnitude))
    
    arr = np.loadtxt(filename)
    
    for i in range(len(magnitude)):
        sig = arr[:,0]
        logM = arr[:,i+1]
        tck = splrep(sig, logM)
        logMmin_new[i] = splev(sigma[i], tck)
            
    return logMmin_new
