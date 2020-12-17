Sequential fitting of HOD parameters using a fast method based on mass binned halo paircounts.

Before running this you must create the halo mass-distance paircount tables in the paircounting folder ../paircounting/

Then change the run_path variable in config_fitting.py to the base path for your paircount tables such that your central-central
paircounts can be found at run_path + "_cencen.npy"

The other variables which will need to be changed in config_fitting are:

###cosmo_factor###

This is the location of the cosmology rescaling factor text file which has been used to correct for the difference between
the Abacus and MXXL cosmologies

###target_2pcf###

This is the location of the text file for the target correlation functions (from MXXL)

###target_num_den###

This contains the target number densities (from MXXL)

###included_params###

This determines which parameters are used in the fitting and which have fixed values taken from the "..._fit.txt" files
This should be set to [1,1,1,1,1] at the start of sequential fitting (So all the parameters are fitted initially)
A 1 means the parameter is being fitted. A 0 means the fixed value is being used.
The order of parameters is [Mmin, sigma_logM, M0, M1, alpha]. This order is used throughout.

The batch script "batch_sequential.sh" contains the pipeline for completing the sequential fitting.
The order of parameter fitting is determined by this script.

###Running the Fitting###

The basic premise is:
1. Fit the HOD parameters to produce the target clustering and number density
2. Fit smooth curves as a function of magnitude to these parameters (Using Alex's code in hod_params)
3. Fix the values of one of the parameters to the smoothed curve
Repeat steps 1.-3. with the fixed parameter no longer included in the fitting

Removing a parameter from the fitting is done by changing the config file using sed. This is sensitive to line number so be careful!
Updating the text files with the smoothed parameter is done by get_fitted_param.py with a command line argument
corresponding to the parameter to be fitted (1 = Mmin, 2 = sigma_logM, etc...)

Altering the order of fitting you simply have to change the sed command to remove parameters in a different order. And also the argument passed
to get_fitted_param.py at each step must be changed to be the parameter which is being removed from the fitting at each step

The batch script should run straight through and on the final step the values of Mmin are adjusted using Alex's code for correcting the number density

This will provide the final parameter values at number_density/params_fixed_number_density.txt

Running check_HOD.py will show the HODs for the final parameters
And running check_2pcf.py will show the final correlation functions and parameters

Let me know if you have any questions - Cameron 17/12/20

