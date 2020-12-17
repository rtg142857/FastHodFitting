import numpy as np
import matplotlib.pyplot as plt


# function for Mmin and M1
def M_function(log_mass, L_s, M_t, a_m):
    # same functional form as Eq. 11 from Zehavi 2011
    # returns magnitude as a function of mass
    lum = L_s*(10**log_mass/M_t)**a_m * np.exp(1-(M_t/10**log_mass))
    magnitudes = 4.76 - 2.5*np.log10(lum)
    return magnitudes


def M0_function(magnitude, A, B, C, D):
    #return 10**(A*magnitude + B)
    return 10**(A*magnitude**3 + B *magnitude**2 + C *magnitude + D)

def sigma_function(magnitude, A, B, C, D):
    return A + (B-A) / (1.+np.exp((magnitude+C)*D))
    #return A*magnitude**3 + B *magnitude**2 + C *magnitude + D


def alpha_function(magnitude, A, B, C):
    return np.log10((A*((4.76-magnitude)/2.5))**B + C)
    #return A*magnitude**3 + B *magnitude**2 + C *magnitude + D



if __name__ == "__main__":

    # plot fits to SDSS HOD parameters
    
    log_mass = np.arange(11,16,0.01)
    plt.plot(M_function(log_mass, 3.9184e+09, 3.0665e+11, 2.5763e-01), log_mass,
             label="Mmin", c="C0")
    plt.plot(M_function(log_mass, 3.7056e+09, 4.7768e+12, 3.0594e-01), log_mass,
             label="M1", c="C3")
    plt.xlim(-17.5, -22.5)
    plt.ylim(11, 16)
    plt.legend(loc="upper left").draw_frame(False)
    plt.xlabel("mag")
    plt.ylabel("M (Msun/h)")
    plt.show()
    
    mags = np.arange(-23,-17,0.1)
    plt.plot(mags, np.log10(M0_function(mags, -0.7134, -2.5844)), label="M0")
    plt.xlim(-17.5, -22.5)
    plt.ylim(9, 14)
    plt.legend(loc="upper left").draw_frame(False)
    plt.xlabel("mag")
    plt.ylabel("M (Msun/h)")
    plt.show()
    
    plt.plot(mags, sigma_function(mags, 0.025836, 0.68127, 21.05, 2.5), label="sigma_logM")
    plt.xlim(-17.5, -22.5)
    plt.ylim(0,0.7)
    plt.legend(loc="upper left").draw_frame(False)
    plt.xlabel("mag")
    plt.ylabel("sigma")
    plt.show()

    plt.plot(mags, alpha_function(mags, 0.098284, 80.2760, 10.0), label="alpha")
    plt.xlim(-17.5, -22.5)
    plt.ylim(0.8, 2.8)
    plt.legend(loc="upper left").draw_frame(False)
    plt.xlabel("mag")
    plt.ylabel("alpha")
    plt.show()
