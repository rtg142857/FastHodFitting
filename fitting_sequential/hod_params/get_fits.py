import numpy as np
import matplotlib.pyplot as plt
import degeneracy as deg
import parameter_fits as fits


def read_params(filename):
    mags, logMmin, sigmalogM, logM0, logM1, alpha = np.loadtxt(filename, unpack=True)
    return mags, logMmin, sigmalogM, logM0, logM1, alpha


def save_params(filename, mags, logMmin, sigmalogM, logM0, logM1, alpha):
    np.savetxt(filename, np.array([mags, logMmin, sigmalogM, logM0, logM1, alpha]).transpose())
    
    
def get_fits(param_file, degeneracy_file, save_file, make_plots=True):
    
    # read best fit HOD parameters
    mags, logMmin, sigma, logM0, logM1, alpha = read_params(param_file)
    
    # Smooth fit to Mmin
    print("Fitting Mmin")
    logMmin_smooth = fits.fit_M(mags, logMmin, p0=[5e9, 3e11, 0.3])

    if make_plots:
        plt.scatter(mags, logMmin, c="C0", edgecolor="none", label="Original fit")
        plt.plot(mags, logMmin_smooth, c="k", ls=":", label="Smoothed")
        plt.title("Mmin")
        plt.legend(loc="upper right").draw_frame(False)
        plt.xlabel("mag")
        plt.ylabel("log(Mmin)")
        plt.show()
        
    # Smooth fit to M1
    print("Fitting M1")
    logM1_smooth = fits.fit_M(mags, logM1, p0=[5e9, 3e11, 0.3])
    
    if make_plots:
        plt.scatter(mags, logM1, c="C0", edgecolor="none", label="Original fit")
        plt.plot(mags, logM1_smooth, c="k", ls=":", label="Smoothed")
        plt.title("M1")
        plt.legend(loc="upper right").draw_frame(False)
        plt.xlabel("mag")
        plt.ylabel("log(M1)")
        plt.show()
        
    # Smooth fit to M0
    print("Fitting M0")
    logM0_smooth = fits.fit_M0(mags, logM0, p0=[-0.007134, 0.05844,-1,1])
    
    if make_plots:
        plt.scatter(mags, logM0, c="C0", edgecolor="none", label="Original fit")
        plt.plot(mags, logM0_smooth, c="k", ls=":", label="Smoothed")
        plt.title("M0")
        plt.legend(loc="upper right").draw_frame(False)
        plt.xlabel("mag")
        plt.ylabel("log(M0)")
        plt.show()
        
    # Smooth fit to alpha
    print("Fitting alpha")
    mask = np.ones(len(mags), dtype="bool")
    #mask[0] = False # ignore the -22 sample
    alpha_smooth = fits.fit_alpha(mags, alpha, mask=mask)
    
    if make_plots:
        plt.scatter(mags, alpha, c="C0", edgecolor="none", label="Original fit")
        plt.plot(mags, alpha_smooth, c="k", ls=":", label="Smoothed")
        plt.title("alpha")
        plt.legend(loc="upper right").draw_frame(False)
        plt.xlabel("mag")
        plt.ylabel("alpha")
        plt.show()
        
    # Transform the new, smoothed values of Mmin to sigma, using the degeneracy
    print("Transforming smoothed Mmin to sigma")
    sigma_trans = deg.logMmin_to_sigma(mags, logMmin_smooth, degeneracy_file)
    
    # Smooth fit to sigma (after transformation)
    print("Fitting sigma")
    mask = np.ones(len(mags), dtype="bool")
    #mask[0], mask[-1], mask[-2] = False, False, False # ignore brightest/faintest samples
    sigma_trans_smooth = fits.fit_sigma(mags, sigma_trans, mask=mask, p0=[0.03, 0.7, 20, 2.5])
    
    if make_plots:
        plt.scatter(mags, sigma, c="C0", edgecolor="none", label="Original fit")
        plt.scatter(mags, sigma_trans, c="C1", edgecolor="none", label="From Mmin")
        plt.plot(mags, sigma_trans_smooth, c="k", ls=":", label="Smoothed")
        plt.title("sigma")
        plt.legend(loc="upper right").draw_frame(False)
        plt.xlabel("mag")
        plt.ylabel("sigma")
        plt.show()
        
    # Transform smoothed values of sigma back to Mmin, using the degeneracy
    logMmin_trans = deg.sigma_to_logMmin(mags, sigma_trans_smooth, degeneracy_file)
    
    # Smooth fit to Mmin (after transformation)
    print("Fitting Mmin")
    logMmin_trans_smooth = fits.fit_M(mags, logMmin_trans, p0=[5e9, 3e11, 0.3])
    
    if make_plots:
        plt.scatter(mags, logMmin, c="C0", edgecolor="none", label="Original fit")
        plt.scatter(mags, logMmin_trans, c="C1", edgecolor="none", label="From sigma")
        plt.plot(mags, logMmin_trans_smooth, c="k", ls=":", label="Smoothed")
        plt.title("Mmin")
        plt.legend(loc="upper right").draw_frame(False)
        plt.xlabel("mag")
        plt.ylabel("Mmin")
        plt.show()
        
    # save new smoothed HOD parameters
    save_params(save_file, mags, logMmin_trans_smooth, sigma_trans_smooth, logM0_smooth, logM1_smooth, alpha_smooth)
    

if __name__ == "__main__":
   
    param_file = "params/fit_params.txt"
    save_file = "params/params_new.txt"
    degeneracy_file = "params/NFW_cut_pg_Sigma_logm_Mmin_fit.txt"

    get_fits(param_file, degeneracy_file, save_file, make_plots=False)
    
