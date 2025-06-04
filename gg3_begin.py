# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 15:35:09 2025

@author: callu
"""

"""
first, consider tracking a single neutron
"""

import random as rnd
import numpy as np

def main():
    
    # constants 
    
    nu = 3 # neutrons liberated per fission
    #mfp = 0.036 # mean free path: distance between interactions
    sigma_f = 1.235 # fission cross section in bn
    sigma_s = 4.566 # scattering cross section in bn
    p_f = sigma_f/(sigma_f + sigma_s) # fission probability
    
    R = 3.6e-11 # Fissioncore sample radius of 36000fm
    #R = 8.37e-2 # real-life sample radius of 8.37cm
    
    d = 2.4e-13 # Fissioncore lattice constant of 240fm
    #d = 2.75e-10 # real-life lattice constant of 275000fm
    
    radius = ((sigma_f + sigma_s)*1e-28/np.pi)**0.5 # U-235 nuclear radius
    
    print("\nnuclear radius:", radius)
    print("\nd - nuclear radius:", d-radius)
    
    collision_nuclei = np.abs(np.sin(np.pi*radius/d)) # collision detection parameter
    
    
    # initialise a neutron with position and velocity i.e x,y,z & v_x,v_y,v_z 
    
    live_neutrons = 1 # number of live neutrons
    v = (2*(2e6)*(1.602e-19)/1.675e-27)**0.5 # neutron speed
    n = [0, 0, 0, 0, 0, 0] # our neutron as an array, initialised at the origin
    
    
    # initialise time and time steps
    
    t = 0 # initialise time
    
    #dt = mfp/v # mean time between interactions
    #dt = 2*radius/v # time for neutron to cross a nucleas
    #dt = 5e-22 # Fissioncore's timestep of 0.5zs - may be too small?
    
    disp = radius # how far we want the neutron to travel each timestep
    dt = disp/v # calculate the associated timestep
    print("\ndt", dt)
    
    steps = int(1e5) # how long will we run for?
    dsteps = 0 # step counter for analysis
    
    
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
    
    position() # update position to prevent collision with the nucleas at the origin

    
    # the main for loop:
    
    for i in range(0, steps):
        # update time and counter
        t += dt
        dsteps += 1
        
        # update position
        position()
        r = (n[0]**2 + n[1]**2 + n[2]**2)**0.5 # neutron distance from sample centre (origin)
        
        # check neutron position to see if a collision has occurred
        x_close = (np.abs(n[0]/d)%1)*d # mod method
        y_close = (np.abs(n[1]/d)%1)*d
        z_close = (np.abs(n[2]/d)%1)*d
        
        """collision_x = np.abs(np.sin(np.pi*n[0]/d)) # sine method
        collision_y = np.abs(np.sin(np.pi*n[1]/d))
        collision_z = np.abs(np.sin(np.pi*n[2]/d))
        
        if ((collision_x <= collision_nuclei) and (collision_y <= collision_nuclei) and (collision_z <= collision_nuclei)):
            velocity() # new velocity vector
            print("COLLISION DETECTED")
            print("dist from origin:", r)
        else:
            print("MISS")"""
        
        # collision detection condition
        if ((x_close <= radius)or(x_close >= d-radius)) and ((y_close <= radius)or(y_close >= d-radius)) and ((z_close <= radius)or(z_close >= d-radius)):
            # if condition above met, neutron is within a cube of a nucleas
            
            nuclei_distance = (x_close**2 + y_close**2 + z_close**2)**0.5 # how close the neutron is to a nuclei
            if nuclei_distance <= radius: # if this distance is less than nuclei radius
            
                velocity() # new velocity vector
            
                prob = rnd.random() # rnd float between 0 and 1
                if prob < p_f: # if fission occurs, create more neutrons
                    live_neutrons += nu - 1 # -1 since bombarding neutron used up in interaction
                    print("FISSION!")
                else: # if scattering occurs, nothing for now
                    print("SCATTER")
            
            #print("\nx_close:", x_close, "\ny_close:", y_close, "\nz_close:", z_close)
                print("Steps:", dsteps)
        
        
        #print("\ncollision_x:", collision_x, "\ncollision_y:", collision_y, "\ncollision_z:", collision_z)
        
    print("\nTime elapsed:", t, "\nLive neutrons:", live_neutrons, "\nDistance from sample centre:", r)

main()

""" 
Review:
    - 
    - 
    - 
"""

#%%
"""
Make a 2-d lattice of uranium-235 nuclei of sample size R
"""
R = 3.6e-11 # sample size of 36000fm
d = 2.4e-13 # lattice constant of 240fm
eta = int(2*R/d - 1) # number of nuclei across sample diameter
print(eta)

x = np.linspace(-2, 2, 5)
y = np.linspace(-2, 2, 5)
print(x)

xmesh, ymesh = np.meshgrid(x,y, sparse=True) # sparse is for efficiency
print(xmesh,ymesh)
coords = []

for (i,j) in zip(xmesh,ymesh):
    for (k,l) in zip(i,j):
        element = (k,l)
        coords.append(element)

print(coords) # all coordinates

#%%
"""
Rough Work
"""
print(int(1e6))

