# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 16:33:24 2025

@author: callu
"""

"""
running the main while loop for 10 trials with the end goal of automating the process of testing
values of sample_radius to obtain a computed value for critical radius
"""

import math
import numpy as np
import random as rnd
import matplotlib.pyplot as plt

from matplotlib.lines import Line2D

def main():

    # material specific values
    
    sigma_f = 1.235 # fission cross section in bn
    sigma_c = 0.089 # capture cross section in bn
    sigma_s = 4.566 + 1.804 # scattering (elastic + inelastic) cross section in bn
    sigma_t = sigma_f + sigma_c + sigma_s # transport cross section
    
    density = 18710 # density of U-235 in kg m^-3
    number_density = 4.794e28 # number density of U-235 in m^-3
    mfp = 1/(sigma_t*1e-28*number_density) # transport mean free path in m
    
    
    sample_radius = 7.43e-2 # real-life sample radius (critical value is 8.37cm)
    
    
    sample_mass = (4*np.pi*density*(sample_radius)**3)/3 # calculation of sample mass
    
    neutron_mass = 1.675e-27 # neutron mass in kg
    k_energy = 2e6 # neutron kinetic energy in eV
    speed = (2*k_energy*(1.602e-19)/neutron_mass)**0.5 # neutron speed in ms^-1
    
    
    # neutrons liberated per fission
    
    nu_mean = 2.637 # average neutrons liberated per fission
    nu_sd = 1 # standard deviation of nu_mean
    nu = 0 # initialise nu value which will be randomly obtained every fission
    nu_0 = 100 # number of initial neutrons
    
    m = 7 # max no. of neutrons emitted (8 if using nu_empprobs)
    
    nu_discrete = [i for i in range(0, m+1)] # storage of integer values nu can take
    nu_empprobs = [0.027, 0.158, 0.339, 0.305, 0.133, 0.038, -0.001, 0.001, 0] # empirical probabilities from MULTIPLICITIES paper
    nu_binomprobs = [] # probabilities from binomial distribution from MULTIPLICITIES paper
    
    for i in range(0, m+1):
        nu_binomprobs.append(math.factorial(m)/(math.factorial(i)*math.factorial(m - i))*((nu_mean/m)**i)*(1-nu_mean/m)**(m-i))


    # initialisation
        
    prob_f = sigma_f/sigma_t # fission probability
    prob_s = sigma_s/sigma_t # scatter probability
    
    t = 0 # initialise running time variable
    palette_choice = 0 # initialise colour choice
        
    
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
        path_length = -mfp*np.log(rnd.random())
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
    
    delta_radius = 0.01 # initialise the change in sample radius
    crit_status = [] # intialise list of status of criticality for each trial
    
    
    try:
        
        for i in range(10): # run the following 10 times (10 trials of simulation)
        
            print("TRIAL {0:}".format(i+1))
            
            # initialise variables and lists for this trial
            
            t = 0 # initialise time for this trial
            live_neutrons = nu_0 # running counter of live neutrons
            escaped_neutrons = 0 # running counter of escaped neutrons
            fisscaptured_neutrons = 0 # running counter of fission captured neutrons
            radcaptured_neutrons = 0 # running counter of radiatively captured neutrons
            scattered_neutrons = 0 # running counter of scattered neutrons

            neutron_box = [[0,0,0,0,0,0,0,0,0] for j in range(nu_0)] # 'box' of empty neutrons: x,y,z, v_x,v_y,v_z, dt, tsort, r
            graph_data = [[t], [live_neutrons]] # list for plotting N(t)
            nu_data = [] # list for plotting nu values
            kappa = [] # list for averaging effective multiplication factor
            
            for neutron in neutron_box:
                velocity(neutron) # initialise velocities
                timestep(neutron) # initialise timesteps
                sortingtime(neutron) # initialise sorting times

        
            # loop to pick neutron with smallest dt and run accordingly
            
            loop = 0 # to keep track of the number of loops
            
            """while live_neutrons > 0 and live_neutrons < 3000 and t < 1e-6: # run until all neutrons escape or limit achieved or sufficient time elapsed"""
            while loop >= 0: # run indefinitely (until conditions below met)
                        
            
                N_m = live_neutrons # number of live neutrons at mth generation
                
                tsort = [neutron[7] for neutron in neutron_box] # a list of the sorting times with same order as neutrons list
                mini = tsort.index((min(tsort))) # the index of the neutron with the smallest sorting time
                
                position(neutron_box[mini]) # update the position of the neutron with the smallest sorting time
                
                t = min(tsort) # the minimum sorting time is the time since initialisation (elapsed time)
                
                if neutron_box[mini][8] > sample_radius: # if neutron escapes sample
                    
                    neutron_box.pop(mini) # remove the specific neutron from the neutrons list
                    live_neutrons -= 1 # account for the escaped neutron
                    escaped_neutrons += 1 # account for the escaped neutron
                     
                else: # if neutron remains inside sample, an interaction will occur
                
                    prob = rnd.random() # rnd float between 0 and 1 to determine fission or scattering
                    
                    if prob <= prob_s: # if scattering occurs, neutron prepares to move off in a new direction
                    
                        velocity(neutron_box[mini]) # give the neutron a new velocity
                        timestep(neutron_box[mini]) # give the neutron a new timestep
                        sortingtime(neutron_box[mini]) # update the neutron's sorting time
                        
                        scattered_neutrons += 1 # neutron scattered
                    
                    elif prob_s < prob <= (prob_s + prob_f): # if fission occurs, capture neutron and create 3 more at collision site
                        
                        # get an integer value for nu this fission
                        #nu = abs(int(round(rnd.gauss(nu_mean, nu_sd)))) # GAUSSIAN distribution with sd = 1, no -ve values
                        #nu = np.random.poisson(nu_mean) # nu value obtained using POISSON statistics
                        basket = rnd.choices(nu_discrete, nu_binomprobs) # basket via pdf based on empirical data from MULTIPLICITIES paper
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
                        
            
                N_mplus1 = live_neutrons # number of live neutrons at m+1th generation
                kappa.append(N_mplus1/N_m) # calculate and add kappa value to list
                
                loop += 1 # loop finished
                
                if loop%10 == 0: # every 10 loops:
                    graph_data[0].append(t) # add this loop's elapsed time value to the data list
                    graph_data[1].append(live_neutrons) # # add this loop's live_neutrons value to the data list
                    
                if loop%10000 == 0: # every 10000 loops:
                    print("\nloop:", loop, "live_neutrons:", live_neutrons) # display no. of live neutrons
        
        
                # conditions to determine the criticality of the trial and break the loop
                
                if live_neutrons == 0:
                    crit_status.append("subcritical") # fizzled out
                    palette_choice = 0 # cyan plots
                    break
                    
                elif live_neutrons > 3000:
                    crit_status.append("supercritical") # assumed indefinite increase
                    palette_choice = 1 # red plots
                    break
    
                elif t > 1e-6: # [not really correct]
                    crit_status.append("critical") # assumed rate constant
                    palette_choice = 2 # olive plots
                    break
    
    
            print(crit_status)
            
    
            # plot data
            
            plt.rcParams['figure.figsize'] = [10, 6] # set the size for plotting
            plt.rcParams['figure.dpi'] = 100 # set the resolution for plotting
            
            palette = ["cyan", "red", "olive"] # 3 colours depending on criticality
            paint = palette[palette_choice] # paint for the painting
            
            plt.plot(graph_data[0], graph_data[1], 
                     color=paint, marker='o', markersize=1) # plot live_neutrons against elapsed time
                            

    except KeyboardInterrupt: # if program terminated, data will still be plotted/printed (code below runs)
        pass

    # plot all the trials together

    plt.axhline(y=3000, color="orange", linestyle="dashed", label="Supercritical Limit") # supercritical condition for live_neutrons

    plt.title("10 trials of live_neutrons vs. t")
    plt.xlabel("t")
    plt.ylabel("live_neutrons")
    plt.figtext(0, 0.045, "sample_radius: {0:.6} cm".format(sample_radius*100))
    plt.figtext(0, 0, "nu_mean: {0:} neutrons".format(nu_mean))
    
    custom_legend = [Line2D([0], [0], color="cyan"), # define the legend
                     Line2D([0], [0], color="red"),
                     Line2D([0], [0], color="olive"),
                     Line2D([0], [0], color="orange", linestyle="dashed")]    
    
    plt.legend(custom_legend, ["Subcritical", "Supercritical", "Critical", "Neutron Limit"], loc="lower right") # legend
    plt.show()
    
    for i in range(10):
        print("TRIAL {0:}: ".format(i+1), crit_status[i]) # show the results of the 10 trials PROBLEM HERE?
        
main()

#%%


treats = [i for i in range(10)]
print(treats)
treats.append(600)
print(treats)
treats = [i for i in range(10)]
print(treats)











