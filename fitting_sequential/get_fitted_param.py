import numpy as np
import sys

# Take the fitted parameter values and fix in a new file to be used in later fitting

labels = ['Magnitude','Mmin','sigma_logM','M0','M1','alpha']

params_new = np.genfromtxt("hod_params/params/params_new.txt")

param_number = int(sys.argv[1])

params = params_new[:,param_number]

np.savetxt(labels[param_number]+'_fit.txt',params)

