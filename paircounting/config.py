import numpy as np

path = "/cosma7/data/durham/dc-grov1/Halo_mass_pair_binning/BGS/reduced_catalog/haloes_z0.2_nfw_cut.hdf5"

boxsize = 2000

r_bin_edges = np.logspace(-2,2.2,85)

mass_bin_edges = np.logspace(11,16,101)
mass_bin_centres = 10**(np.log10(mass_bin_edges[:-1])+(np.diff(np.log10(mass_bin_edges))/2))
num_sat_parts = 3

run_label = "subsample_pairs_test"

def var_sample(a,b): # slope = a, full sample mass = b
    subsampling = 10**(a*(np.log10(mass_bin_centres) -  b))
    subsampling[subsampling > 1] = 1
    return subsampling

print(mass_bin_edges)
print(mass_bin_centres)

subsample_array = var_sample(2,12)

