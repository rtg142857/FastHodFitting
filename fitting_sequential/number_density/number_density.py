import numpy as np
from hod_5param import HOD_5param


def match_target_density(log_Mmin, sigma, log_M0, log_M1, alpha, log_n_targ, log_mass_halo, log_n_halo):
    """
    Returns the value of log_Mmin needed to match a target number density
    Args:
        log_Mmin:      log10 of halo mass with 50% probability of hosting a central galaxy [Msun/h]
        sigma:         Width of step in central HOD
        log_M0:        log10 of halo mass where there is a cutoff in the satellite HOD [Msun/h]
        log_M1:        log10 of average halo mass that hosts 1 satellite galaxy [Msun/h]
        alpha:         Power law slope for satellite HOD
        log_n_targ:    log10 of target number density [Mpc/h]^{-3}
        log_mass_halo: array of log10 of halo mass bins [Msun/h]
        log_n_halo:    array of log10 of halo number densities [Mpc/h]^{-3}
    """
    
    log_Mmin_min, log_Mmin_max, log_Mmin_mid = log_Mmin-1, log_Mmin+1, log_Mmin
    
    hod = HOD_5param(10**log_Mmin_mid, sigma, 10**log_M0, 10**log_M1, alpha)
    log_n_mid = np.log10(hod.number_density(log_mass_halo, log_n_halo))
    
    while abs(log_n_mid - log_n_targ) > 1e-5:
        
        log_n_mid_prev = log_n_mid
        
        if log_n_targ < log_n_mid:
            log_Mmin_min = log_Mmin_mid
        else:
            log_Mmin_max = log_Mmin_mid
            
        log_Mmin_mid = (log_Mmin_min + log_Mmin_max) / 2.

        hod = HOD_5param(10**log_Mmin_mid, sigma, 10**log_M0, 10**log_M1, alpha)
        log_n_mid = np.log10(hod.number_density(log_mass_halo, log_n_halo))
        
    return log_Mmin_mid



if __name__ == "__main__":
    
    # HOD parameters from fits with sigma_logM, M1 and M0 fixed, for -22 sample
    log_Mmin = 1.366664735943966669e+01 
    sigma    = 6.145474399999999449e-01 
    log_M0   = 1.369773217406736876e+01 
    log_M1   = 1.443356999983215516e+01 
    alpha    = 1.298748031241174328e+00
    
    # read file of mass function measured from Abacus simulation
    log_mass_halo, log_n_halo = np.loadtxt("mass_function.dat", unpack=True)

    # read file with target number densities
    magnitudes, log_n_targs = np.loadtxt("target_number_density.dat", unpack=True)
    
    # get Mmin needed to match target number density
    log_Mmin_new = match_target_density(log_Mmin, sigma, log_M0, log_M1, alpha, log_n_targs[0], log_mass_halo, log_n_halo)
    
    print("Original log_Mmin:\t", log_Mmin)
    print("Match number density:\t", log_Mmin_new)
    print("Ratio:\t\t\t", 10**log_Mmin_new/10**log_Mmin)
