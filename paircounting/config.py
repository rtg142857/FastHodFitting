import numpy as np

path = "/cosma7/data/durham/dc-grov1/Halo_mass_pair_binning/BGS/haloes_z0.2.hdf5"

boxsize = 2000

r_bin_edges = np.logspace(-2,2.2,85)

mass_bin_edges = np.logspace(1,6,101)

num_sat_parts = 3

run_label = "test_0.1"

subsample_array = np.ones(100)*0.1 # Constant 10% for now
