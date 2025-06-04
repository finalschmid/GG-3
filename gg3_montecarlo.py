# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 17:28:36 2025

@author: callu
"""

"""
Using a randomly sampled path length instead of a lattice of nuclei
"""

import random as rnd
import numpy as np

def main():
    
    mfp = 0.036 # mean free path
    nu = 3 # neutrons liberated per fission
    R = 8.37e-2 # real-life sample radius of 8.37cm
    
    sigma_f = 1.235 # fission cross section in bn
    sigma_s = 4.566 # scattering cross section in bn
    p_f = sigma_f/(sigma_f + sigma_s) # fission probability
    
    live_neutrons = 1 # number of live neutrons
    v = (2*(2e6)*(1.602e-19)/1.675e-27)**0.5 # neutron speed REMAINS CONSTANT
    n = [0, 0, 0, 0, 0, 0] # our neutron as an array, initialised at the origin

    steps = int(1e3) # how long will we run for?
    dsteps = 0 # step counter for analysis

    t = 0 # initialise time
    path_length = 0 # how far we want the neutron to travel each timestep VARIES
    dt = 0 # initialise timestep VARIES
    """path_length = -mfp*np.log(rnd.random()) # putting this here to see if it needs to be initialised like this 
    dt = path_length/v"""
    
    
    # velocity function: gives the neutron a random velocity direction
    
    def velocity():
        theta = rnd.uniform(0, 2*np.pi) # random theta value on sphere
        phi = rnd.uniform(0, np.pi) # random phi value on sphere
        
        n[3] = v*np.cos(theta)*np.sin(phi) # v_x
        n[4] = v*np.sin(theta)*np.sin(phi) # v_y
        n[5] = v*np.cos(phi) # v_z
        
    velocity() # call for a random initial direction to move in
    
    
    # position function: updates the neutron's position
    
    def position():
        n[0] += n[3]*dt # update x-position 
        n[1] += n[4]*dt # update y-position
        n[2] += n[5]*dt # update z-position
        
        
    # the main for loop:
    
    for i in range(0, steps):
        # decide on path length and timestep
        path_length = -mfp*np.log(rnd.random()) # path length from exponential probability density
        dt = path_length/v # calculate the associated timestep
        
        # update time and counter
        t += dt
        dsteps += 1
        
        # update position
        position()

        # calculate distance from sample centre
        r = (n[0]**2 + n[1]**2 + n[2]**2)**0.5
        
        # check if the neutron has escaped the sample
        if r > R: # if the neutron outside the sample
            live_neutrons -= 1 # account for the escaped neutron
            print("**NEUTRON ESCAPED!**")
            break # stop the for loop once the neutron has escaped

        # decide on the type of interation
        prob = rnd.random() # rnd float between 0 and 1
        if prob < p_f: # if fission occurs, create more neutrons
            live_neutrons += nu - 1 # -1 since bombarding neutron used up in interaction
            print("FISSION!")
        else: # if scattering occurs, nothing for now
            print("SCATTER")

        # next direction
        velocity()
        
        # check on our neutron
        #print("\nn:", n)
        
    print("\nSteps Elapsed:", dsteps,"\nTime elapsed:", t, "\nLive neutrons:", live_neutrons, "\nDistance from sample centre:", r)
    
main()
#%%

a = 2
a -= 1

print(a)

