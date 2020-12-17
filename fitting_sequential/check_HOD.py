import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import config
import config_fitting
# Load parameters from config file
mass_bin_edges = config.mass_bin_edges

# Load parameters from fitting config file

Cen_HOD = config_fitting.Cen_HOD
spline_kernel_integral = config_fitting.spline_kernel_integral
cumulative_spline_kernel = config_fitting.cumulative_spline_kernel
Sat_HOD = config_fitting.Sat_HOD

if __name__=="__main__":
    """
    Small script to check the HODs are being implemented correctly
    Take fitted HOD parameters and check the plot looks
    as expected
    Change the parameter array in config_fitting.py to [1,1,1,1,1] before running this script
    """
    plt.figure(figsize=(8,8))
    params = np.genfromtxt("number_density/params_fixed_number_density.txt")[:,1:]
    Mmins = params[:,0]
    sigma_logMs = params[:,1]
    M0s = params[:,2]
    M1s = params[:,3]
    alphas = params[:,4]
    for i in range(9):
        param = np.genfromtxt("number_density/params_fixed_number_density.txt")[i,1:]
        print(param)
        cen = Cen_HOD(param,mass_bin_edges,Mmins[i],sigma_logMs[i])
        sat = Sat_HOD(param,cen,mass_bin_edges,M0s[i],M1s[i],alphas[i])
        sat[np.isnan(sat)] = 0
    #plt.plot(mass_bin_edges,cen,label = "Central")
    #plt.plot(mass_bin_edges,sat,label = "Satellite")
        plt.plot(mass_bin_edges,cen+sat,label = "Total")
    plt.ylim(bottom = 1e-4) # Low number cutoff
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Halo Mass /Solar Masses")
    plt.ylabel("Number of galaxies per halo")
    plt.legend()
    plt.show()
    plt.close()
    


