# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 16:00:05 2025

@author: callu
"""

"""
for plotting miscellaneous data of interest
"""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

def main():

    # choose between U-235 (0) & Pu-239 (1)
    
    ert = 0


    # material specific values for U-235 & Pu-239
    
    sigma_f = [1.235, 1.8] # fission cross section in bn
    sigma_c = [0.089, 0.053] # capture cross section in bn
    sigma_s = [4.566 + 1.804, 4.394 + 1.46] # scattering (elastic + inelastic) cross section in bn
    sigma_t = [sigma_f[i] + sigma_c[i] + sigma_s[i] for i in range(2)]# transport cross section
    
    number_density = [4.794e28, 3.93e28] # number densities in m^-3
    mfp = 1/(sigma_t[ert]*1e-28*number_density[ert]) # transport mean free path in m


    # data for plotting
    
    r = np.arange(0.001, 1, 0.001) # excluding 0 as log would diverge and 1 as rnd method does not include it
    l = -mfp*np.log(r)*100 # calculate paths in cm


    # plotting
        
    plt.rcParams['figure.figsize'] = [5, 3] # set the size for plotting
    plt.rcParams['figure.dpi'] = 100 # set the resolution for plotting

    plt.axvline(x=mfp*100, color="orange", linestyle="dashed", label="$\lambda_t = {0:.2}$cm".format(mfp*100)) # mfp plotted
    #plt.plot(r, l, color="purple", linestyle="-", label="$\ell$") # equation plotted
    plt.plot(l, r, color="purple", linestyle="-", label="$\ell$") # FLIP AXES

    
    plt.title("Neutron Path Length Distribution: $\ell = -\lambda_t \ln r$")
    plt.ylabel(r"Random Number r")
    #plt.xlim(0,1)
    plt.xlim(0,) # FLIP AXES
    plt.xlabel("Path Length $\ell$ (cm)")
    #plt.ylim(0,) # thanks cian :)
    plt.ylim(0,1) # FLIP AXES
    
    plt.annotate("$\infty$", xytext=(15, 0.2), xy=(19, 0.2),
            arrowprops=dict(arrowstyle="->"), fontsize=15)
    #plt.figtext(0.8, 0.375, "to infinity")
    
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.show()

main()

#%%
ert = 1

sigma_f = [1.235, 1.8] # fission cross section in bn
sigma_c = [0.089, 0.053] # capture cross section in bn
sigma_s = [4.566 + 1.804, 4.394 + 1.46] # scattering (elastic + inelastic) cross section in bn
sigma_t = [sigma_f[i] + sigma_c[i] + sigma_s[i] for i in range(2)] # transport cross section


prob_f = sigma_f[ert]/sigma_t[ert] # fission probability
prob_s = sigma_s[ert]/sigma_t[ert] # scatter probability

# neutron speed calulation

neutron_mass = 1.675e-27 # neutron mass in kg
k_energy = 2e6 # neutron kinetic energy in eV
speed = (2*k_energy*(1.602e-19)/neutron_mass)**0.5 # neutron speed in ms^-1

#print("{0:.2e}".format(speed))
#density = [18710, 15600] # mass densities in kg m^-3
density = [18923, 19800] # densities used in papers

core_radius = 4.9538e-2 # core radius (critical values are 8.37cm and 6.346cm)
core_mass = (4*np.pi*density[ert]*(core_radius)**3)/3 # calculation of core mass

print("{0:}".format(core_mass))
#%%

"""
Plotting the nu histograms
"""
import math
import numpy as np
import random as rnd
import matplotlib.pyplot as plt

ert = 0
material_core = ["Uranium-235", "Plutonium-239"] # core materials

# neutrons liberated per fission

nu_mean = [2.637, 3.172] # average neutrons liberated per fission
nu_sd = 1 # standard deviation of nu_mean
nu = 0 # initialise nu value which will be randomly obtained every fission
nu_0 = 100 # number of initial neutrons

nu_data = [] # list for plotting nu values
nu_neg = [] # list for seeing how many negative values there are for the gaussian case

m = [5, 6] # max no. of neutrons emitted (8 if using nu_empprobs)

nu_discrete = [i for i in range(0, m[ert]+1)] # storage of integer values nu can take
nu_empprobs = [0.027, 0.158, 0.339, 0.305, 0.133, 0.038, -0.001, 0.001, 0] # empirical probabilities from MULTIPLICITIES paper
nu_binomprobs = [] # probabilities from binomial distribution from MULTIPLICITIES paper

for i in range(0, m[ert]+1):
    nu_binomprobs.append(math.factorial(m[ert])/(math.factorial(i)*math.factorial(m[ert] - i))*((nu_mean[ert]/m[ert])**i)*(1-nu_mean[ert]/m[ert])**(m[ert]-i))

for i in range(10000):
    # get an integer value for nu this fission
    """nu = int(round(rnd.gauss(nu_mean[ert], nu_sd))) # GAUSSIAN distribution with sd = 1, no -ve values"""
    """nu = np.random.poisson(nu_mean[ert]) # nu value obtained using POISSON statistics"""
    basket = rnd.choices(nu_discrete, nu_binomprobs) # basket via pdf based on EMPIRICAL/BINOMIAL data from MULTIPLICITIES paper
    nu = basket[0] # take the integer from the basket"""
    nu_data.append(nu) # add the above nu value to the list
    
    if nu < 0: # keeping track of negative values from Gaussian case
        nu_neg.append(nu)


# plot histogram of nu values

plt.rcParams['figure.figsize'] = [10, 6] # set the size for plotting THE HISTOGRAM
plt.rcParams['figure.dpi'] = 100 # set the resolution for plotting
plt.rcParams.update({'font.size': 15}) # set font size

plt.hist(nu_data, color="lightgreen", ec="black", bins=len(set(nu_data))) # number of distinct nu values = number of bins"""

plt.axvline(x=nu_mean[ert], color="black", linestyle="dashed", linewidth=3,
            label=r"$\bar\nu ={0:}$".format(nu_mean[ert])) # supercritical condition for live_neutrons"""
"""plt.axvline(x=2.47, color="red", linestyle="dashed",linewidth=3,
label=r"$\bar\nu = 2.47$") # EMPIRICAL"""

plt.title(r"Neutrons liberated per fission $\nu$ (BINOMIAL)") # r in front of string to make a "raw" string
plt.xlabel(r"$\nu$ value") # renders \nu instead of \n u
plt.ylabel("frequency")

"""plt.figtext(0, 0.035, "core_radius: {0:.3} cm".format(core_radius*100))
plt.figtext(0, 0, "tamp_radius: {0:} cm".format(tamp_radius*100))"""
plt.figtext(0.7, 0.035, "core: {0:}".format(material_core[ert]), fontsize=20)
"""plt.figtext(0.7, 0, "tamper: {0:}".format(material_tamp[0]))"""

plt.legend(fontsize=20)
plt.tight_layout()
plt.show()

print("\nOccurrences of -1 neutrons: {0:} times".format(len(nu_neg)),
      "\nbins:", len(set(nu_data)),
      "\nAvg: {0:}".format(sum(nu_data)/len(nu_data)))

