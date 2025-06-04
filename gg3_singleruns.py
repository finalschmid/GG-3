# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 18:09:37 2025

@author: callu
"""

"""
tracking multiple neutrons and using a monte carlo method to randomly sample path length as suggested 
in "An Introduction to Computer Simulation Methods" by Gould & Tobochnik.

Used for single simulation runs and also includes the tamper code.
Use tamped_untamped.py for tamped runs.
"""

import math
import numpy as np
import random as rnd
import matplotlib.pyplot as plt

#from IPython.display import set_matplotlib_formats
#set_matplotlib_formats('svg')

def main():
    
    # choose between U-235 (0) & Pu-239 (1)
    
    ert = 1
    material_core = ["Uranium-235", "Plutonium-239"] # core materials
    

    # material specific values for U-235 & Pu-239
    
    sigma_f = [1.235, 1.8] # fission cross section in bn
    sigma_c = [0.089, 0.053] # capture cross section in bn
    sigma_s = [4.566 + 1.804, 4.394 + 1.46] # scattering (elastic + inelastic) cross section in bn
    sigma_t = [sigma_f[i] + sigma_c[i] + sigma_s[i] for i in range(2)] # transport cross section
    
    mass_number = [235.04, 239.05] # Uranium and Plutonium mass numbers
    density = [18710, 15600] # mass densities in kg m^-3
    number_density = [4.794e28, 3.93e28] # number densities in m^-3
    mfp_core = 1/(sigma_t[ert]*1e-28*number_density[ert]) # core transport mean free path in m
    
    
    # material specific values for the tamper (DU for now)
    
    mfp_tamp = 4.342e-2 # (Reed tb) value for tamper mfp in m
    density_tamp = 18950 # (Reed tb) mass density in kg m^-3
    
    
    # sample radii and masses
    
    core_radius = 5.55e-2 # core radius (critical values are 8.37cm and 6.346cm)
    
    tamp_radius = 0 # tamper radius, i.e. shell thickness
    core_mass = (4*np.pi*density[ert]*(core_radius)**3)/3 # calculation of core mass
    tamp_mass = (4*np.pi*density_tamp*((core_radius + tamp_radius)**3 - (core_radius)**3))/3 # calculation of tamper mass
    
    
    # neutron speed calulation
    
    neutron_mass = 1.675e-27 # neutron mass in kg
    k_energy = 2e6 # neutron kinetic energy in eV
    speed = (2*k_energy*(1.602e-19)/neutron_mass)**0.5 # neutron speed in ms^-1
    
    
    # neutrons liberated per fission
    
    nu_mean = [2.637, 3.172] # average neutrons liberated per fission
    """nu_sd = 1 # standard deviation of nu_mean"""
    nu = 0 # initialise nu value which will be randomly obtained every fission
    nu_0 = 100 # number of initial neutrons
    
    m = [5, 6] # max no. of neutrons emitted (8 if using nu_empprobs)
    
    nu_discrete = [i for i in range(0, m[ert]+1)] # storage of integer values nu can take
    #nu_empprobs = [0.027, 0.158, 0.339, 0.305, 0.133, 0.038, -0.001, 0.001, 0] # empirical probabilities from MULTIPLICITIES paper
    nu_binomprobs = [] # probabilities from binomial distribution from MULTIPLICITIES paper
    
    for i in range(0, m[ert]+1):
        nu_binomprobs.append(math.factorial(m[ert])/(math.factorial(i)*math.factorial(m[ert] - i))*((nu_mean[ert]/m[ert])**i)*(1-nu_mean[ert]/m[ert])**(m[ert]-i))


    # calculation of time for core consumption

    mfp_f = 1/(sigma_f[ert]*1e-28*number_density[ert]) # fission mean free path
    nuclei = 1e3*core_mass*6.02e23/mass_number[ert] # nuclei in core
    generations = np.log(nuclei*(nu_mean[ert]-1)/nu_0)/np.log(nu_mean[ert]) # generations of neutron liberation

    t_fiss = generations*mfp_f/speed # time for all fissile material to be consumed
    #t_fiss = t_fiss*0.5 # 1/2
    

    # initialisation
    
    live_neutrons = nu_0 # running counter of live neutrons
    total_neutrons = nu_0 # running counter of total neutrons that have been in the simulation
    interactions = 0 # running counter of interactions which have been simulated
    
    escaped_neutrons = 0 # running counter of escaped neutrons
    fisscaptured_neutrons = 0 # running counter of fission captured neutrons
    radcaptured_neutrons = 0 # running counter of radiatively captured neutrons
    scattered_neutrons = 0 # running counter of scattered neutrons
    
    prob_f = sigma_f[ert]/sigma_t[ert] # fission probability
    prob_s = sigma_s[ert]/sigma_t[ert] # scatter probability
    
    t = 0 # initialise running time
    mfp = mfp_core # initialise live mfp as mfp_core because neutrons start at core centre
    palette_choice = 0 # initialise colour choice
    
    neutron_box = [[0,0,0,0,0,0,0,0,0] for j in range(nu_0)] # 'box' of empty neutrons: x,y,z, v_x,v_y,v_z, dt, tsort, r
    graph_data = [[t], [live_neutrons]] # list for plotting N(t)
    nu_data = [] # list for plotting nu values
    
    
    # function to give random velocity
    
    def velocity(i):
        # i refers to the specific neutron
        theta = rnd.uniform(0, 2*np.pi) # random theta value on sphere
        phi = rnd.uniform(0, np.pi) # random phi value on sphere
            
        i[3] = speed*np.cos(theta)*np.sin(phi) 
        i[4] = speed*np.sin(theta)*np.sin(phi) 
        i[5] = speed*np.cos(phi)
    
    
    # function to give random timestep to a neutron based on mfp
    
    def timestep(i):
        # i refers to the specific neutron
        path_length = -mfp*np.log(rnd.random()) # conditions will change this mfp value
        dt = path_length/speed
            
        i[6] = dt # timestep, used in position calculation
    
    
    # function to update the sorting time of a neutron
    
    def sortingtime(i):
        i[7] += i[6] # tsort, used to find the 'earliest' neutron
    
       
    # function to update the position and calculate distance of a neutron
    
    def position(i):
        # i refers to the specific neutron
        i[0] += i[3]*i[6] # x = v_x*dt
        i[1] += i[4]*i[6] # y = v_y*dt
        i[2] += i[5]*i[6] # z = v_z*dt
    
        i[8] = (i[0]**2 + i[1]**2 + i[2]**2)**0.5 # r = sqrt(x^2 + y^2 + z^2)
        
    
    # initialise neutrons
    
    for neutron in neutron_box:
        velocity(neutron) # initialise velocities
        timestep(neutron) # initialise timesteps
        sortingtime(neutron) # initialise sorting times
        
    
    # loop to pick neutron with smallest dt and run accordingly
    
    loop = 0 # to keep track of the number of loops
    crit_status = 0 # intialise the status of criticality
    
    try:
        
        while loop >= 0: # run indefinitely (until conditions below met)
              
            """N_m = live_neutrons # number of live neutrons at mth generation"""
            
            tsort = [neutron[7] for neutron in neutron_box] # a list of the sorting times with same order as neutrons list
            mini = tsort.index((min(tsort))) # the index of the neutron with the smallest sorting time
            
            position(neutron_box[mini]) # update the position of the neutron with the smallest sorting time
            
            t = min(tsort) # the minimum sorting time is the time since initialisation (elapsed time)
            
            if neutron_box[mini][8] > (core_radius + tamp_radius): # if neutron escapes sample
                
                neutron_box.pop(mini) # remove the specific neutron from the neutrons list
                live_neutrons -= 1 # account for the escaped neutron
                escaped_neutrons += 1 # account for the escaped neutron
                 
            elif neutron_box[mini][8] > core_radius: # if neutron inside tamper, a scatter will occur
            
                mfp = mfp_tamp # mfp inside tamper
                
                velocity(neutron_box[mini]) # give the neutron a new velocity
                timestep(neutron_box[mini]) # give the neutron a new timestep
                sortingtime(neutron_box[mini]) # update the neutron's sorting time
                
                scattered_neutrons += 1 # neutron scattered

            else: # otherwise, the neutron is in the core
                
                interactions += 1 # an interaction occurs
                mfp = mfp_core # mfp inside core
                prob = rnd.random() # rnd float between 0 and 1 to determine fission or scattering
                
                if prob <= prob_s: # if scattering occurs, neutron prepares to move off in a new direction
                
                    velocity(neutron_box[mini]) # give the neutron a new velocity
                    timestep(neutron_box[mini]) # give the neutron a new timestep
                    sortingtime(neutron_box[mini]) # update the neutron's sorting time
                    
                    scattered_neutrons += 1 # neutron scattered
                
                elif prob_s < prob <= (prob_s + prob_f): # if fission occurs, capture neutron and create 3 more at collision site
                    
                    # get an integer value for nu this fission
                    """nu = abs(int(round(rnd.gauss(nu_mean[ert], nu_sd)))) # GAUSSIAN distribution with sd = 1, no -ve values"""
                    """nu = np.random.poisson(nu_mean) # nu value obtained using POISSON statistics"""
                    basket = rnd.choices(nu_discrete, nu_binomprobs) # basket via pdf based on EMPIRICAL/BINOMIAL data from MULTIPLICITIES paper
                    nu = basket[0] # take the integer from the basket
                    nu_data.append(nu) # add the above nu value to the list
                    
                    # fission site coordinates (distance unecessary as position function will be called again later)
                    fiss_x = neutron_box[mini][0] # these are used in the following for loop
                    fiss_y = neutron_box[mini][1]
                    fiss_z = neutron_box[mini][2]
                    
                    # for loop to add nu more neutrons at fission site with new velocities and timesteps
                    for i in range(nu): # do this for each new neutron (i.e nu times)
                        neutron_box.append([fiss_x, fiss_y, fiss_z, 0, 0, 0, 0, min(tsort), 0]) # add to end of neutron_box
                        
                        velocity(neutron_box[-1]) # run for last neutron in neutron_box
                        timestep(neutron_box[-1]) 
                        sortingtime(neutron_box[-1])
                        
                    neutron_box.pop(mini) # bombarding neutron is at same index in neutron_box because of append
                    
                    live_neutrons += nu - 1 # -1 since bombarding neutron used up in fission
                    total_neutrons += nu 
                    fisscaptured_neutrons += 1 # bombarding neutron used up in fission
                    
                else: # if capture occurs, neutron removed from list
                    
                    neutron_box.pop(mini) # remove the captured neutron
                    
                    live_neutrons -= 1 # 1 neutron captured
                    radcaptured_neutrons += 1 # update counter
                    
        
            """N_mplus1 = live_neutrons # number of live neutrons at m+1th generation
            kappa.append(N_mplus1/N_m) # calculate and add kappa value to list"""
            
            loop += 1 # loop finished
            
            if loop%10 == 0: # every 10 loops:
                graph_data[0].append(t) # add this loop's elapsed time value to the data list
                graph_data[1].append(live_neutrons) # # add this loop's live_neutrons value to the data list
                
            if loop%10000 == 0: # every 10000 loops:
                print("\nTime Elapsed: {0:.3f} μs".format(t*1e6), # display time and no. of live neutrons
                      "live_neutrons:", live_neutrons) 
    
    
            # conditions to determine criticality and break the loop
            
            if live_neutrons == 0: # if neutrons run out
                crit_status = "subcritical" # fizzled out
                palette_choice = 0 # cyan plots
                break
                
            elif t > t_fiss: # if sample used up
                crit_status = "supercritical" # assumed indefinite increase
                palette_choice = 1 # red plots
                break

    
    except KeyboardInterrupt: # if program terminated, data will still be plotted/printed (code below runs)
        crit_status = "supercritical" # putting this here for termination of long runs
        palette_choice = 1 
        pass


    # set parameters globally
    
    plt.rcParams['figure.figsize'] = [10, 6] # set the size for plotting THE HISTOGRAM
    plt.rcParams['figure.dpi'] = 300 # set the resolution for plotting
    plt.rcParams.update({'font.size': 10}) # set font size
    
    plt.rcParams['xtick.labelsize'] = 15 # size of axis numbers
    plt.rcParams['ytick.labelsize'] = 15
    
    palette = ["cyan", "red", "olive"] # 3 colours depending on criticality
    paint = palette[palette_choice] # paint for the painting
    palette_ln = ["blue", "firebrick", "olive"]
    paint_ln = palette_ln[palette_choice] # paint for the ln(painting)


    # plot live_neutrons against elapsed time
    
    plt.plot(graph_data[0], graph_data[1], color=paint, marker='o', markersize=1, 
             label="$N$ Time: {0:.3f} μs".format(graph_data[0][-1]*1e6))
    plt.axvline(x=t_fiss, color="orange", linestyle="dashed", linewidth=2,
                label="Time Limit: {0:.3f} μs".format(t_fiss*1e6)) # supercritical condition for live_neutrons
    
    plt.title("Live Neutrons $N$ vs. Elapsed Time t", fontsize=15)
    plt.xlabel("t (s)", fontsize=15)
    plt.ylabel("Live Neutrons $N$", fontsize=15)
    
    plt.xlim(0,)
    plt.ylim(0,)
    
    plt.figtext(0, 0.035, "Sample Radius: {0:.3} cm".format(core_radius*100), fontsize=15)
    plt.figtext(0.65, 0.035, "Sample: {0:}".format(material_core[ert]), fontsize=15)
    
    plt.figtext(0.8, 0.95, "(LINEAR)", fontsize=15)

    plt.legend(fontsize=15)
    plt.tight_layout()
    plt.grid()
    plt.show()
    
    
    # plot the natural logarithm of live_neutrons against elapsed time
    
    plt.plot(graph_data[0], np.log(graph_data[1]), color=paint_ln, marker="o", markersize=1, 
             label="$\ln(N)$ Time: {0:.3f} μs".format(graph_data[0][-1]*1e6))
    plt.axvline(x=t_fiss, color="orange", linestyle="dashed", linewidth=2,
                label="Time Limit: {0:.3f} μs".format(t_fiss*1e6)) # supercritical condition for live_neutrons

    """linefit = np.poly1d(np.polyfit(graph_data[0], np.log(graph_data[1]), 1)) # straight line equation fitted to the data
    ax2.plot(graph_data[0], linefit(graph_data[0]), "g", label="Fitted Line") # plot the straight line"""
    
    plt.title("$\ln$(Live Neutrons) $\ln(N)$ vs. Elapsed Time t", fontsize=15)
    plt.xlabel("t (s)", fontsize=15)
    plt.ylabel("$\ln$(Live Neutrons) $\ln(N)$", fontsize=15)
    
    plt.xlim(0,)
    plt.ylim(0,)
    
    plt.figtext(0, 0.035, "Sample Radius: {0:.3} cm".format(core_radius*100), fontsize=15)
    plt.figtext(0.65, 0.035, "Sample: {0:}".format(material_core[ert]), fontsize=15)
    
    plt.figtext(0.8, 0.95, "(LOGARITHMIC)", fontsize=15)

    plt.legend(fontsize=15)
    plt.grid()
    plt.tight_layout() # apply "tight_layout" to figure
    plt.show()
    
    
    # plot histogram of nu values
    
    """nu_avg = (sum(nu_data))/(len(nu_data))
    
    plt.hist(nu_data, color="lightgreen", ec="black", bins=len(set(nu_data))) # number of distinct nu values = number of bins 
    #plt.axvline(x=nu_mean[ert], color="black", linestyle="dashed", linewidth=3,
            #label=r"$\bar\nu ={0:}$".format(nu_mean[ert])) # supercritical condition for live_neutrons
    plt.axvline(x=nu_avg, color="red", linestyle="dashed", linewidth=3,
            label=r"$\bar\nu ={0:}$".format(nu_mean[ert])) # supercritical condition for live_neutrons

    plt.title(r"Neutrons liberated per fission $\nu$", fontsize=20) # r in front of string to make a "raw" string
    plt.xlabel(r"$\nu$ value", fontsize=20) # renders \nu instead of \n u
    plt.ylabel("frequency", fontsize=20)
    
    plt.ylim(0,)
    
    plt.figtext(0, 0.035, "Sample Radius: {0:.3} cm".format(core_radius*100), fontsize=15)
    plt.figtext(0.65, 0.035, "Sample: {0:}".format(material_core[ert]), fontsize=15)
    
    plt.grid()
    plt.legend(fontsize=20)
    plt.tight_layout()
    plt.show()"""
    
    
    # plot pie chart of interactions
    
    pie_labels = ["scatter", "fission", "capture"]
    sigma_labels = ["$\sigma_{scatter}$", "$\sigma_{fission}$", "$\sigma_{capture}$"]
    
    pie_sectors = [scattered_neutrons, fisscaptured_neutrons, radcaptured_neutrons]
    sigma_sectors = [sigma_s[ert], sigma_f[ert], sigma_c[ert]]
    #hatches = ["x", "o", "."]
    
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.set_aspect("equal")
    ax1.set_aspect("equal")

    ax1.pie(pie_sectors, labels=pie_labels, autopct="%1.1f%%", textprops={'fontsize': 15})
    ax1.set_title(r"Simulated Interactions" "\n" r"({0:})".format(material_core[ert]), fontsize=15)

    ax2.pie(sigma_sectors, labels=sigma_labels, autopct="%1.1f%%", textprops={'fontsize': 15},
            colors=["palegreen", "lightcoral", "gold"])
    ax2.set_title("{0:} cross-sections".format(material_core[ert]), fontsize=15)

    plt.show()
    
    # print data
    
    """alpha = linefit.c[0] # 1d polynomial 1st order coefficient (i.e the slope)
    kappa_bar = sum(kappa)/len(kappa) # mean kappa value"""
    
    
    print(
        "\n",
        "\nCore Mass: {0:.2f} kg".format(core_mass),
        "\nCore Radius: {0:.2f} cm".format(core_radius*100),
        "\nTamper Mass: {0:.2f} kg".format(tamp_mass),
        "\nTamper Radius: {0:.2f} cm".format(tamp_radius*100),
        "\n",
        "\nTime for Core Consumption: {0:.3f} μs".format(t_fiss*1e6),
        "\nTime Elapsed: {0:.3f} μs".format(t*1e6),
        "\nCriticality: {0:}".format(crit_status),
        "\n",
        "\nInteractions Simulated:", interactions,
        "\nNeutrons Simulated:", total_neutrons, 
        "\nNeutrons Escaped:", escaped_neutrons,
        "\n",
        "\nFission Events:", fisscaptured_neutrons,
        "\nRadiative Capture Events:", radcaptured_neutrons,
        "\nScattering Events:", scattered_neutrons,
          )
    

main()

#%%


ert = 1

sigma_f = [1.235, 1.8] # fission cross section in bn
sigma_c = [0.089, 0.053] # capture cross section in bn
sigma_s = [4.566 + 1.804, 4.394 + 1.46] # scattering (elastic + inelastic) cross section in bn
sigma_t = [sigma_f[i] + sigma_c[i] + sigma_s[i] for i in range(2)] # transport cross section

prob_f = sigma_f[ert]/sigma_t[ert] # fission probability
prob_s = sigma_s[ert]/sigma_t[ert] # scatter probability
prob_c = sigma_c[ert]/sigma_t[ert] # capture probability

print(prob_s, prob_f, prob_c, 
      prob_s + prob_f + prob_c)


pie_labels = ["scatter", "fission", "capture"]
pie_sectors = [sigma_s[ert], sigma_f[ert], sigma_c[ert]]

plt.pie(pie_sectors, labels=pie_labels, autopct="%1.1f%%")
plt.title("Likelihood of Interactions ", fontsize=15)

plt.tight_layout()
plt.show()




