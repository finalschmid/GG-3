# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 20:56:58 2025

@author: callu
"""

"""
for working on the collision detection, keeping values simple for debugging
"""

import random as rnd
import numpy as np

t = 0 #time
dt = 1 #timestep
steps = 100 # number of timesteps

d = 0.5 # lattice constant
radius = 0.05 # nuclei radius

v = 0.25 # neutron speed
n = [0, 0, 0, 0, 0, 0] # our neutron

def velocity(): # give random velocity direction
        theta = rnd.uniform(0, 2*np.pi) # random theta value on sphere
        phi = rnd.uniform(0, np.pi) # random phi value on sphere
        
        n[3] = v*np.cos(theta)*np.sin(phi) # v_x
        n[4] = v*np.sin(theta)*np.sin(phi) # v_y
        n[5] = v*np.cos(phi) # v_z

velocity()

def position(): # update position
    n[0] += n[3]*dt # update x-position 
    n[1] += n[4]*dt # update y-position
    n[2] += n[5]*dt # update z-position
    
for i in range(0, steps):
    t += dt
    position()
    r = (n[0]**2 + n[1]**2 + n[2]**2)**0.5 # distance from origin
    
    #print("\nneutron:", n, "\ndistance:", r)
    
    x_close = (np.abs(n[0]/d)%1)*d
    y_close = (np.abs(n[1]/d)%1)*d
    z_close = (np.abs(n[2]/d)%1)*d
    
    #print("\nx_close:", x_close, "\ny_close:", y_close, "\nz_close:", z_close)
    
    if ((x_close <= radius)or(x_close >= d-radius)) and ((y_close <= radius)or(y_close >= d-radius)) and ((z_close <= radius)or(z_close >= d-radius)):
        velocity() # new velocity vector
        print("COLLISION DETECTED")
        print("\nx_close:", x_close, "\ny_close:", y_close, "\nz_close:", z_close)
        print("\ntime:", t)
    
    
    
#%%

def rng01():
    prob = rnd.random() # rnd float between 0 and 1
    print(prob)
    
for i in range(0,5):
    rng01()
    
    
    
    
    