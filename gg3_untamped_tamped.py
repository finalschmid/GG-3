# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 21:52:09 2025

@author: callu
"""

"""
Use this for comparing an untamped and then tamped core
"""

import math
import numpy as np
import random as rnd
import matplotlib.pyplot as plt

from matplotlib.lines import Line2D

def main():

    # choose between U-235 (0) & Pu-239 (1)
    
    ert = 0 # choose core
    #tre = 1 # choose tamper
    material_core = ["Uranium-235", "Plutonium-239"] # core materials
    material_tamp = ["DU", "BeO"] # tamper materials

    # material specific values for U-235 & Pu-239
    
    sigma_f = [1.235, 1.8] # fission cross section in bn
    sigma_c = [0.089, 0.053] # capture cross section in bn
    sigma_s = [4.566 + 1.804, 4.394 + 1.46] # scattering (elastic + inelastic) cross section in bn
    sigma_t = [sigma_f[i] + sigma_c[i] + sigma_s[i] for i in range(2)]# transport cross section
    
    mass_number = [235.04, 239.05] # Uranium and Plutonium mass numbers
    density = [18710, 15600] # mass densities in kg m^-3
    number_density = [4.794e28, 3.93e28] # number densities in m^-3
    mfp_core = 1/(sigma_t[ert]*1e-28*number_density[ert]) # core transport mean free path in m


    # material specific values for the tamper 
    
    mfp_tamp = [4.342e-2, 2.541e-2] # (Reed tb) value for tamper mfp in m
    density_tamp = [18950, 3020] # (Reed tb) mass density in kg m^-3
    
    
    # sample radii and masses
    
    core_radius = 7e-2 # core radius (critical values are 8.37cm and 6.346cm)
    tamp_radius = 1e-2 # tamper radius, i.e. shell thickness
    
    core_mass = (4*np.pi*density[ert]*(core_radius)**3)/3 # calculation of core mass
    tamp_mass1 = (4*np.pi*density_tamp[0]*((core_radius + tamp_radius)**3 - (core_radius)**3))/3 # calculation of tamper mass
    tamp_mass2 = (4*np.pi*density_tamp[1]*((core_radius + tamp_radius)**3 - (core_radius)**3))/3 
    
    # neutron speed calulation
    
    neutron_mass = 1.675e-27 # neutron mass in kg
    k_energy = 2e6 # neutron kinetic energy in eV
    speed = (2*k_energy*(1.602e-19)/neutron_mass)**0.5 # neutron speed in ms^-1
    
    
    # neutrons liberated per fission
    
    nu_mean = [2.637, 3.172] # average neutrons liberated per fission
    nu = 0 # initialise nu value which will be randomly obtained every fission
    nu_0 = 100 # number of initial neutrons
    
    graph_data_t = [] # list where each run's list of t values will be stored
    graph_data_live_neutrons = []

    m = [5, 6] # max no. of neutrons emitted (8 if using nu_empprobs)

    nu_discrete = [i for i in range(0, m[ert]+1)] # storage of integer values nu can take
    nu_binomprobs = [] # probabilities from binomial distribution from MULTIPLICITIES paper
    
    for i in range(0, m[ert]+1):
        nu_binomprobs.append(math.factorial(m[ert])/(math.factorial(i)*math.factorial(m[ert] - i))*((nu_mean[ert]/m[ert])**i)*(1-nu_mean[ert]/m[ert])**(m[ert]-i))


    # calculation of time for core consumption

    mfp_f = 1/(sigma_f[ert]*1e-28*number_density[ert]) # fission mean free path
    nuclei = 1e3*core_mass*6.02e23/mass_number[ert] # nuclei in core
    generations = np.log(nuclei*(nu_mean[ert]-1)/nu_0)/np.log(nu_mean[ert]) # generations of neutron liberation

    t_fiss = generations*mfp_f/speed # time for all fissile material to be consumed
    #t_fiss = t_fiss/2 # HALF

    # initialisation
        
    prob_f = sigma_f[ert]/sigma_t[ert] # fission probability
    prob_s = sigma_s[ert]/sigma_t[ert] # scatter probability
    
    t = 0 # initialise running time variable
    mfp = mfp_core # initialise live mfp as mfp_core because neutrons start at core centre
    #palette_choice = 0 # initialise colour choice
        
    
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
        path_length = -mfp*np.log(rnd.random()) # random method is half open range 0 <= r < 1
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
    
    
    # automation parameters
    
    trials = 3 # untamped then tamped
    tamp_radius_live = 0 # intitialise tamper radius used in loop
    crit_status = [] # intialise list of status of criticality for each trial
    term_t = [] # initialise list of termination times for the 2 runs
    #paint_board = [] # list of paints for corresponding trial
    #paint_board_ln = []

    
    try:
        
        for i in range(trials): # run the following trials times
        
            this_trial = i # for keeping track of which trial this is later
            if this_trial > 0: # if this is after the first trial
                tamp_radius_live = tamp_radius # run the tamped version
                tamper_trial = i - 1
                print(mfp_tamp[tamper_trial])
            print("TRIAL {0:}".format(i+1))
            
            
            # initialise variables and lists for this trial
            
            t = 0 # initialise time for this trial
            live_neutrons = nu_0 # running counter of live neutrons
            escaped_neutrons = 0 # running counter of escaped neutrons
            fisscaptured_neutrons = 0 # running counter of fission captured neutrons
            radcaptured_neutrons = 0 # running counter of radiatively captured neutrons
            scattered_neutrons = 0 # running counter of scattered neutrons

            neutron_box = [[0,0,0,0,0,0,0,0,0] for j in range(nu_0)] # 'box' of empty neutrons: x,y,z, v_x,v_y,v_z, dt, tsort, r
            graph_data_t.append([t]) # add this loop's list for t values
            graph_data_live_neutrons.append([live_neutrons])
            nu_data = [] # list for plotting nu values
            
            for neutron in neutron_box:
                velocity(neutron) # initialise velocities
                timestep(neutron) # initialise timesteps
                sortingtime(neutron) # initialise sorting times

        
            # loop to pick neutron with smallest dt and run accordingly
            
            loop = 0 # to keep track of the number of loops
            
            while loop >= 0: # run indefinitely (until conditions below met)
                            
                tsort = [neutron[7] for neutron in neutron_box] # a list of the sorting times with same order as neutrons list
                mini = tsort.index((min(tsort))) # the index of the neutron with the smallest sorting time
                
                position(neutron_box[mini]) # update the position of the neutron with the smallest sorting time
                
                t = min(tsort) # the minimum sorting time is the time since initialisation (elapsed time)
                
                if neutron_box[mini][8] > (core_radius + tamp_radius_live): # if neutron escapes sample
                    
                    neutron_box.pop(mini) # remove the specific neutron from the neutrons list
                    live_neutrons -= 1 # account for the escaped neutron
                    escaped_neutrons += 1 # account for the escaped neutron
                
                elif neutron_box[mini][8] > core_radius: # if neutron inside tamper, a scatter will occur
                    
                    mfp = mfp_tamp[tamper_trial] # mfp inside tamper
                    
                    velocity(neutron_box[mini]) # give the neutron a new velocity
                    timestep(neutron_box[mini]) # give the neutron a new timestep
                    sortingtime(neutron_box[mini]) # update the neutron's sorting time
                    
                    scattered_neutrons += 1 # neutron scattered

                else: # otherwise, the neutron is in the core
                
                    mfp = mfp_core # mfp inside core
                    prob = rnd.random() # rnd float between 0 and 1 to determine fission or scattering
                    
                    if prob <= prob_s: # if scattering occurs, neutron prepares to move off in a new direction
                    
                        velocity(neutron_box[mini]) # give the neutron a new velocity
                        timestep(neutron_box[mini]) # give the neutron a new timestep
                        sortingtime(neutron_box[mini]) # update the neutron's sorting time
                        
                        scattered_neutrons += 1 # neutron scattered
                    
                    elif prob_s < prob <= (prob_s + prob_f): # if fission occurs, capture neutron and create 3 more at collision site
                        
                        # get an integer value for nu this fission
                        """nu = abs(int(round(rnd.gauss(nu_mean, nu_sd)))) # GAUSSIAN distribution with sd = 1, no -ve values
                        nu = np.random.poisson(nu_mean) # nu value obtained using POISSON statistics"""
                        basket = rnd.choices(nu_discrete, nu_binomprobs) # basket via pdf based on EMPIRICAL data/BINOMIAL from MULTIPLICITIES paper
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
                        fisscaptured_neutrons += 1 # bombarding neutron used up in fission
                        
                    else: # if capture occurs, neutron removed from list
                        
                        neutron_box.pop(mini) # remove the captured neutron
                        
                        live_neutrons -= 1 # 1 neutron captured
                        radcaptured_neutrons += 1 # update counter
                        
            
                #N_mplus1 = live_neutrons # number of live neutrons at m+1th generation
                
                loop += 1 # loop finished
                
                if loop%10 == 0: # every 10 loops:
                    graph_data_t[this_trial].append(t) # add this loop's (i) elapsed time value to this loops list
                    graph_data_live_neutrons[this_trial].append(live_neutrons)
                    
                if loop%10000 == 0: # every 10000 loops:
                    print("\nTime Elapsed: {0:.3f} μs".format(t*1e6), # display time and no. of live neutrons
                      "live_neutrons:", live_neutrons)
        
        
                # conditions to determine the criticality of the trial and break the loop
                
                if live_neutrons == 0: # if neutrons run out
                    crit_status.append("subcritical") # fizzled out
                    #palette_choice = 0 # cyan plots
                    
                    term_t.append(graph_data_t[this_trial][-1]) # time of termination
                    break
                    
                elif t > t_fiss: # if sample used up
                    crit_status.append("supercritical") # assumed indefinite increase
                    #palette_choice = 1 # red plots
                    
                    term_t.append(graph_data_t[this_trial][-1]) # time of termination
                    break
        
    
            # store paint
            
            palette = ["cyan", "lime", "lightgray"] # 2 colours for untamped and tamped
            #paint_board.append(palette[palette_choice]) # paint for the paint_board
            palette_ln = ["blue", "forestgreen", "black"]
            #paint_board_ln.append(palette_ln[palette_choice])
                                        
            
    except KeyboardInterrupt: # if program terminated, data will still be plotted/printed (code below runs)
        pass


    # set parameters globally

    plt.rcParams['figure.figsize'] = [10, 6] # set the size for plotting (10,6 default)
    plt.rcParams['figure.dpi'] = 300 # set the resolution for plotting

    plt.rcParams['xtick.labelsize'] = 15 # size of axis numbers
    plt.rcParams['ytick.labelsize'] = 15


    # linear plots

    for i in range(trials):
        plt.plot(graph_data_t[i], graph_data_live_neutrons[i], 
                     color=palette[i], marker='o', markersize=1) # plot live_neutrons vs. elapsed time
        plt.annotate("{0:}".format(i+1), (graph_data_t[i][-1], graph_data_live_neutrons[i][-1]), 
                     textcoords="offset points", xytext=(10,10), ha='center', 
                     arrowprops=dict(arrowstyle="->"))

    plt.axvline(x=t_fiss, color="orange", linestyle="dashed", linewidth=2) # supercritical condition for live_neutrons

    plt.title("Untamped vs. Tamped core", fontsize=15)
    plt.xlabel("t (s)", fontsize=15)
    plt.ylabel("Live Neutrons $N$", fontsize=15)
    
    plt.xlim(0,)
    plt.ylim(0,)
        
    plt.figtext(0, 0, "Core Radius: {0:.3} cm".format(core_radius*100), fontsize=12)
    plt.figtext(0.2, 0, "Tamper Radius: {0:.3} cm".format(tamp_radius*100), fontsize=12)
    plt.figtext(0.6, 0, "Core: {0:}".format(material_core[ert]), fontsize=12)
    plt.figtext(0.8, 0, "Tampers: {0:}, {1:}".format(material_tamp[0], material_tamp[1]), fontsize=12)

    plt.figtext(0.8, 0.95, "(LINEAR)", fontsize=15)

    
    custom_legend = [Line2D([0], [0], color="cyan"), # define the legend
                     Line2D([0], [0], color="lime"),
                     Line2D([0], [0], color="lightgray"),
                     Line2D([0], [0], color="orange", linestyle="dashed")]    
    plt.legend(custom_legend, 
               ["Untamped Time: {0:.3f} μs".format(term_t[0]*1e6), 
                "DU Time: {0:.3f} μs".format(term_t[1]*1e6),
                "BeO Time: {0:.3f} μs".format(term_t[2]*1e6), 
                "Time Limit: {0:.3f} μs".format(t_fiss*1e6)], 
                fontsize=15) # legend
    plt.grid()
    plt.tight_layout()
    plt.show()
    
    
    # remove diverging logs (replace ln(0) with ln(1))
    
    for i in range(trials): # for each trial
    
        if graph_data_live_neutrons[i][-1] == 0: # if the last element is 0
            graph_data_live_neutrons[i].pop(-1) # remove it
            graph_data_live_neutrons[i].append(1) # replace it with 1
                
    
    # log plots
    
    for i in range(trials):
        plt.plot(graph_data_t[i], np.log(graph_data_live_neutrons[i]), 
                     color=palette_ln[i], marker='o', markersize=1)
        plt.annotate("{0:}".format(i+1), (graph_data_t[i][-1], np.log(graph_data_live_neutrons[i][-1])), 
                     textcoords="offset points", xytext=(10,10), ha='center', 
                     arrowprops=dict(arrowstyle="->"))

    plt.axvline(x=t_fiss, color="orange", linestyle="dashed", linewidth=2) 

    plt.title("Untamped vs. Tamped core", fontsize=15)
    plt.xlabel("t (s)", fontsize=15)
    plt.ylabel("$\ln$(Live Neutrons) $\ln(N)$)", fontsize=15)
    
    plt.xlim(0,)
    plt.ylim(0,)
    
    plt.figtext(0, 0, "Core Radius: {0:.3} cm".format(core_radius*100), fontsize=12)
    plt.figtext(0.2, 0, "Tamper Radius: {0:.3} cm".format(tamp_radius*100), fontsize=12)
    plt.figtext(0.6, 0, "Core: {0:}".format(material_core[ert]), fontsize=12)
    plt.figtext(0.8, 0, "Tampers: {0:}, {1:}".format(material_tamp[0], material_tamp[1]), fontsize=12)
    
    plt.figtext(0.8, 0.95, "(LOGARITHMIC)", fontsize=15)
    
    custom_legend = [Line2D([0], [0], color="blue"), 
                     Line2D([0], [0], color="forestgreen"),
                     Line2D([0], [0], color="black"),
                     Line2D([0], [0], color="orange", linestyle="dashed")]    
    
    plt.legend(custom_legend, 
               ["Untamped Time: {0:.3f} μs".format(term_t[0]*1e6), 
                "DU TIme: {0:.3f} μs".format(term_t[1]*1e6), 
                "BeO Time: {0:.3f} μs".format(term_t[2]*1e6), 
                "Time Limit: {0:.3f} μs".format(t_fiss*1e6)], 
               fontsize=15)
    plt.grid()
    plt.tight_layout()
    plt.show()


    # printing
            
    print(
          "\nCore Mass: {0:.2f} kg".format(core_mass),
          "\nCore Radius: {0:.2f} cm".format(core_radius*100),
          "\nDU Mass: {0:.2f} kg".format(tamp_mass1),
          "\nBeO Mass: {0:.2f} kg".format(tamp_mass2),
          "\nTamper Radius: {0:.2f} cm".format(tamp_radius*100),
          "\nInitiators Used: {0:} neutrons".format(nu_0),
          "\n",
          "\nSUBCRITICAL: {0:} times".format(crit_status.count("subcritical")),
          "\nSUPERCRITICAL: {0:} times".format(crit_status.count("supercritical"))
          )
    
main()