import numpy as np
import config_fitting
save_path = config_fitting.save_path
target_num_den = config_fitting.target_num_den
import matplotlib.pyplot as plt
Cen_HOD = config_fitting.Cen_HOD
Sat_HOD = config_fitting.Sat_HOD

mass_bin_edges = np.logspace(10,16,300) # Originally had 120 bins, but 30 good enough
mass_bin_centres = 10**(np.log10(mass_bin_edges[:-1])+(np.diff(np.log10(mass_bin_edges))/2))

def plot_HODs(HODs):
    """
    Plot HOD/s resulting from fits
    """
    plt.figure(figsize = (8,8))
    for i in range(len(HODs[0,:])):
        plt.plot(mass_bin_edges[:-1] + np.diff(mass_bin_edges)/2,HODs[:,i],c="C"+str(i),label=target_num_den[i,0])
    plt.ylim(1e-3)
    plt.yscale("log")
    plt.xscale("log")
    plt.legend(title = "Magnitude")
    plt.ylabel("n")
    plt.xlabel("Halo Mass /Solar Masses")
    plt.savefig(save_path+"_smooth_HODs.png",bbox_inches="tight")
    plt.close()
    return 0



def get_HODs(params):
    hod_cen_big = Cen_HOD(params,mass_bin_centres)
    hod_sat_big = Sat_HOD(params,hod_cen_big,mass_bin_centres)
    return hod_cen_big + hod_sat_big


HOD_params = np.genfromtxt(save_path+"_best_params.txt")

HODs = get_HODs(HOD_params)

plot_HODs(HODs)