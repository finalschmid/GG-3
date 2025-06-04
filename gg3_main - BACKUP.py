# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 18:09:37 2025

@author: callu
"""

"""
tracking multiple neutrons and using a monte carlo method to randomly sample path length as suggested 
in "An Introduction to Computer Simulation Methods" by Gould & Tobochnik.
"""

import numpy as np
import random as rnd
import matplotlib.pyplot as plt

def main():
    
    mfp = 0.036 # transport mean free path for U-235
    sample_radius = 8e-2 # real-life sample radius (critical value is 8.37cm)
    
    nu = 3 # neutrons liberated per fission
    nu_0 = 100 # number of initial neutrons
    
    live_neutrons = nu_0 # running counter of live neutrons
    escaped_neutrons = 0 # running counter of escaped neutrons
    fissioned_neutrons = 0 # running counter of fissioned neutrons
    
    speed = (2*(2e6)*(1.602e-19)/1.675e-27)**0.5 # neutron speed
    
    t = 0 # initialise running time
    
    sigma_f = 1.235 # fission cross section in bn
    sigma_s = 4.566 # scattering cross section in bn
    prob_fiss = sigma_f/(sigma_f + sigma_s) # fission probability
    
    neutron_box = [[0,0,0,0,0,0,0,0,0] for j in range(nu_0)] # 'box' of empty neutrons: x,y,z, v_x,v_y,v_z, dt, tsort, r
    graph_data = [[t], [live_neutrons]] # list for plotting N(t)
    
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
        
    
    # initialise neutrons
    
    for neutron in neutron_box:
        velocity(neutron) # initialise velocities
        timestep(neutron) # initialise timesteps
        sortingtime(neutron) # initialise sorting times
        
    
    # loop to pick neutron with smallest dt and run accordingly
    
    loop = 0 # to keep track of the number of loops
    
    try:
        
        while live_neutrons > 0 and live_neutrons < 10000: #and loop < 160000: # run until all neutrons escape or limit achieved
                    
            tsort = [neutron[7] for neutron in neutron_box] # a list of the sorting times with same order as neutrons list
            mini = tsort.index((min(tsort))) # the index of the neutron with the smallest sorting time
            
            position(neutron_box[mini]) # update the position of the neutron with the smallest sorting time
            
            t = min(tsort) # the minimum sorting time is the time since initialisation (elapsed time)
            
            if neutron_box[mini][8] > sample_radius: # if neutron escapes sample
                #print("\n**NEUTRON ESCAPED!**", "r:", neutron_box[mini][8])
                
                neutron_box.pop(mini) # remove the specific neutron from the neutrons list
                live_neutrons -= 1 # account for the escaped neutron
                escaped_neutrons += 1 # account for the escaped neutron
                 
            else: # if neutron remains inside sample, an interaction will occur
            
                prob = rnd.random() # rnd float between 0 and 1 to determine fission or scattering
                
                if prob < prob_fiss: # if fission occurs, neutron dies and 3 more will be created at collision site
                    #print("\nFISSION!")
                    
                    # fission site coordinates (distance unecessary as position function will be called again later)
                    fiss_x = neutron_box[mini][0] # these are used in the following for loop
                    fiss_y = neutron_box[mini][1]
                    fiss_z = neutron_box[mini][2]
                    
                    # for loop to add nu more neutrons at fission site with new velocities and timesteps
                    for i in range(nu): # do this for each new neutron (i.e nu times)
                        """neutron_box.insert(i, [fiss_x, fiss_y, fiss_z, 0, 0, 0, 0, min(tsort), 0]) # add nu empty neutrons at indices 0,1,2,...nu-1"""
                        neutron_box.append([fiss_x, fiss_y, fiss_z, 0, 0, 0, 0, min(tsort), 0]) # add to end of neutron_box
                        
                        """velocity(neutron_box[i]) # give the neutrons velocity
                        timestep(neutron_box[i]) # give the neutrons timesteps
                        sortingtime(neutron_box[i]) # update the neutrons' sorting times"""
                        
                        velocity(neutron_box[-1]) # run for last neutron in neutron_box
                        timestep(neutron_box[-1]) 
                        sortingtime(neutron_box[-1])
                        
                    """neutron_box.pop(mini + nu) # remove bombarding neutron (which has increased in index by nu) from neutrons list"""
                    neutron_box.pop(mini) # bombarding neutron is at same index in neutron_box because of append
                    
                    live_neutrons += nu - 1 # -1 since bombarding neutron used up in fission
                    fissioned_neutrons += 1 # bombarding neutron used up in fission
                    
                else: # if scattering occurs, neutron prepares to move off in a new direction
                    #print("\nSCATTER")
                
                    velocity(neutron_box[mini]) # give the neutron a new velocity
                    timestep(neutron_box[mini]) # give the neutron a new timestep
                    sortingtime(neutron_box[mini]) # update the neutron's sorting time
        
            loop += 1 # loop finished
            
            if loop%100 == 0: # every 100 loops:
                graph_data[0].append(t) # add this loop's elapsed time value to the data list
                graph_data[1].append(live_neutrons) # # add this loop's live_neutrons value to the data list
                
            if loop%1000 == 0: # every 1000 loops:
                print("\nloop:", loop, "live_neutrons:", live_neutrons) # display no. of live neutrons
    
    except KeyboardInterrupt: # if program terminated, data will still be plotted
        pass

    # plot data

    plt.plot(graph_data[0], graph_data[1], "ro") # plot live_neutrons against elapsed time
    
    plt.title("live_neutrons vs. t")
    plt.xlabel("t")
    plt.ylabel("live_neutrons")
    plt.figtext(0, 0.045, "sample_radius: {0:.3} cm".format(sample_radius*100))
    plt.figtext(0, 0, "nu: {0:} neutrons".format(nu))
    
    plt.show()
    
    plt.plot(graph_data[0], np.log(graph_data[1]), "bo") # plot the natural logarithm of live_neutrons against elapsed time
    
    linefit = np.poly1d(np.polyfit(graph_data[0], np.log(graph_data[1]), 1)) # straight line equation fitted to the data
    plt.plot(graph_data[0], linefit(graph_data[0]), "g") # plot the straight line
    
    plt.title("ln(live_neutrons) vs. t")
    plt.xlabel("t")
    plt.ylabel("ln(live_neutrons)")
    plt.figtext(0, 0.045, "sample_radius: {0:.3} cm".format(sample_radius*100))
    plt.figtext(0, 0, "nu: {0:} neutrons".format(nu))
    
    plt.show()
    
    
    # print data
    
    alpha = linefit.c[0] # 1d polynomial 1st order coefficient (i.e the slope)
    
    print(
        "\nelapsed time:", t,
        "\n",
        "\nlive_neutrons:", live_neutrons, 
        "\nescaped_neutrons:", escaped_neutrons,
        "\nfissioned_neutrons:", fissioned_neutrons,
        "\n",
        #"\nlinefit:", linefit,
        "\nalpha:", alpha
          )
    
    """for neutron in neutron_box:
        print("\n",neutron) # to see the effects"""

main()

#%%


foo = [
       [0,0,0],
       [0,0,0]
       ]

#foo = [[neutron[parameter] + 2 for parameter in neutron] for neutron in foo]
for neutron in foo: neutron[0] = 2


for i in range(2): foo.insert(i, [0,0,0,0,0,0,0,0,0])
#print(foo)

a = [
     [0,3],
     [8,4],
     [9,1]
     ]

#for i in a: print(min(i[1]))

214721455
74517088




