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
    Input your HOD parameters on the next line and check the plot looks
    as expected
    """
    params = np.array([12.,0.2,10.,13.,1.3 ])
    cen = Cen_HOD(params,mass_bin_edges)
    sat = Sat_HOD(params,cen,mass_bin_edges)

    plt.figure(figsize=(8,8))
    plt.plot(mass_bin_edges,cen,label = "Central")
    plt.plot(mass_bin_edges,sat,label = "Satellite")
    plt.plot(mass_bin_edges,cen+sat,label = "Total")
    plt.ylim(bottom = 1e-4) # Low number cutoff
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Halo Mass /Solar Masses")
    plt.ylabel("Number of galaxies per halo")
    plt.legend()
    plt.show()
    plt.close()
    


