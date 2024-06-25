import numpy as np
import yaml

path_config_filename = sys.argv[1] # Config file path

with open(path_config_filename, "r") as file:
    path_config = yaml.safe_load(file)



path = "../"
path_tracer = "../"
boxsize = 2000

# for wp(rp) only
wp_flag = True # If True calculate wp(rp), otherwise xi(r)
pi_max = 80 # pi max when calculating wp(rp), in Mpc/h
d_pi = 1    # width of pi bins (must divide pi_max)
z_snap = path_config["Params"]["redshift"] # snapshot redshift



# cosmology
flamingo_param_file_path = path_config["Paths"]["params_path"]
with open(flamingo_param_file_path, "r") as file:
    flamingo_params = yaml.safe_load(file)

Om0 = flamingo_params["Cosmology"]["Omega_cdm"]
Ol0 = flamingo_params["Cosmology"]["Omega_lambda"]
#Om0 = 0.3151918
#Ol0 = 0.6848082

r_bin_edges = np.logspace(-2,2,25) # either rp or r, depending if we are calculating wp(rp) or xi(r)

# Due to the units in the file, mass bins here are in units of 10^10 Solar Masses,
# This changes to Solar masses in the fitting so adjust the mass bins there accordingly
mass_bin_edges = np.logspace(0,6,31) # Originally had 120 bins, but 30 good enough
mass_bin_centres = 10**(np.log10(mass_bin_edges[:-1])+(np.diff(np.log10(mass_bin_edges))/2))
num_sat_parts = 3

run_label = "pairs"

def var_sample(a,b): # slope = a, complete sample mass = b
    subsampling = 10**(a*(10 + np.log10(mass_bin_centres) -  b))
    subsampling[subsampling > 1] = 1
    return subsampling


subsample_array = var_sample(2,12)
#subsample_array = np.ones(len(mass_bin_centres))
