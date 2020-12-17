#! /usr/bin/env python
import numpy as np
from scipy.special import erf
import spline

class HOD_5param():

    def __init__(self, Mmin, sigma_logM, M0, M1, alpha, spline=True):
        """
        5 parameter HOD
        Args:
            Mmin: Halo mass with 50% probability of hosting a central galaxy [Msun/h]
            sigma_logM: Width of step in central HOD
            M0:   Halo mass where there is a cutoff in the satellite HOD [Msun/h]
            M1:   Average halo mass that hosts 1 satellite galaxy [Msun/h]
            alpha: Power law slope for satellite HOD
            [spline]: If True, uses a pseudo-Gaussian spline kernel to define the central HOD. 
                      If False, uses a standard Gaussian function. Default value is True. 
        """
        
        self.Mmin = Mmin
        self.M1 = M1
        self.M0 = M0
        self.sigma_logM = sigma_logM
        self.alpha = alpha
        
        self.__spline = spline

    
    def number_centrals_mean(self, log_mass):
        """
        Average number of central galaxies in each halo
        Args:
            log_mass:  array of the log10 of halo mass (Msun/h)
        Returns:
            array of mean number of central galaxies
        """

        if self.__spline:
            # use pseudo Gaussian spline kernel
            return spline.cumulative_spline_kernel(log_mass, 
                    mean=np.log10(self.Mmin), sig=self.sigma_logM/np.sqrt(2))
        else:
            # use error function
            return 0.5*(1 + erf((log_mass - np.log10(self.Mmin))/self.sigma_logM))


    def number_satellites_mean(self, log_mass):
        """
        Average number of satellite galaxies in each halo
        Args:
            log_mass:  array of the log10 of halo mass (Msun/h)
        Returns:
            array of mean number of satellite galaxies
        """
        
        num_cent = self.number_centrals_mean(log_mass)
        num_sat = num_cent * ((10**log_mass - self.M0)/self.M1)**self.alpha
        num_sat[np.where(np.isnan(num_sat))[0]] = 0
        return num_sat


    def number_galaxies_mean(self, log_mass):
        """
        Average total number of galaxies in each halo
        Args:
            log_mass:  array of the log10 of halo mass (Msun/h)
        Returns:
            array of mean number of galaxies
        """
        return self.number_centrals_mean(log_mass) + self.number_satellites_mean(log_mass)

    
    
    def number_density(self, log_mass, log_n_halo):
        """
        Returns number density of galaxies predicted by this HOD, for a given halo mass function
        Args:
            log_mass: array of log10 of halo mass bins (Msun/h)
            n_halo:   array of log10 of halo number densities (Mpc/h)^{-3}
        Returns:
            number density of galaxies in (Mpc/h)^{-3}
        """
        
        N_gal = self.number_galaxies_mean(log_mass)
        n_gal = np.sum(10**log_n_halo * N_gal)
        return n_gal