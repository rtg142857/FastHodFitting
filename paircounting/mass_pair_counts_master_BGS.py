# Import modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import halotools
import halotools.mock_observables
import h5py
import sys
import time
import Corrfunc
from Corrfunc.theory.DD import DD
import fasthod
print("reading in data")

# Load data from specific halo catalogue
path = sys.argv[1]
num_mass_bins = int(sys.argv[2])
run_label = sys.argv[3]

# In this case we have an hdf5 file so read in using h5py


x, y, z, Mvir, is_central, halo_id = fasthod.read_hdf5(path)

x, y, z, Mvir, x_sat, y_sat, z_sat, Mvir_sat = fasthod.split_cen_sat(x,y,z,Mvir,is_central,halo_id)

mass_bin_edges=np.logspace(1,6,num_mass_bins+1)

r_bin_edges = np.logspace(-2,2.2,85)


boxsize = 2000


import multiprocessing

num_threads = 2  # multiprocessing.cpu_count()
print("CPUS total ",num_threads)
print("total particles ",len(x))

start_time = time.time()
samples_test = fasthod.mass_mask(x,y,z,Mvir,mass_bin_edges)
end_time_1 = time.time()
print('starting pair counting')
npairs_test = fasthod.create_npairs_corrfunc(samples_test,samples_test,r_bin_edges,boxsize,num_threads)
npairs_mass_r_bins_test = fasthod.npairs_conversion(samples_test,samples_test,npairs_test,r_bin_edges)
print('pair counting done')
end_time_2 = time.time()
np.save(run_label+".npy",npairs_mass_r_bins_test)

print(run_label)
print(path)
print(num_mass_bins)
print("time to split into mass bins: ",end_time_1 - start_time)
print("total time to run code including paircounting: ",end_time_2 - start_time)
