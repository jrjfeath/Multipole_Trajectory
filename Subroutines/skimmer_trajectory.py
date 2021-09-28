import random 

import numpy as np

from math import sin, cos

def positive_or_negative():
    '''
    Returns a random positive or negative.
    '''
    return 1 if random.random() < 0.5 else -1

def skimmer_trajectory(d):
    '''
    Determine the trajectory of molecules out of the skimmer.\n
    Returns a dictionary containing the velocity in the x, y, and z directions.\n
    Additionally returns the off axis position of the molecules.\n
    velocity = {vx,vy,vz,ry,rz}
    '''
    #Calculate a random gaussian velocity based around theoretical velocity
    grnd = np.random.normal(0, 1) #gaussian random +/- 2.5
    urnd = random.uniform(0.0, 1.0) #uniform random distribution
    grnda = grnd * d['fwhmsk'] #Calculate angle of molecules
    vel = d['velocity'] + (grnd * d['fwhmv']) #Calculate velocity towards collision
    
    #Calculate a uniform float between 0,1 to simulate an even distribution
    #of velocity angles, multiply by 2pi to get -/+ angles
    urnda = urnd * d['tpi']
    #Calculate velocities in the y & z directions
    vy = vel * sin(grnda) * sin(urnda) 
    vz = vel * sin(grnda) * cos(urnda)
    
    ry, rz = d.get('skmr_pos') #Grab off-axis positions
    dist = d.get('skmr_dist') #Check what kind of distribution was selected
    skmr_r = d.get('skmr_radius') * 1000 #Grab skimmer radius

    #Calculate a Gaussian centred at 0 with 90% of the integral in [-sr,sr]
    #We need to divide the radius by 2 to get even proportions in -/+
    if dist == 'g':
        mu, sigma = 0, (skmr_r / 2 / 1.644854)
        ry, rz = np.random.normal(mu, sigma, 2) / 1000

    #Calculate a uniform distribution centred about 0 with limits of skmr_r
    if dist == 'u':
        ry = (random.uniform(0.0, skmr_r) * positive_or_negative()) / 1000
        rz = (random.uniform(0.0, skmr_r) * positive_or_negative()) / 1000

    velocity = {
        "vx" : vel,
        "vy" : vy,
        "vz" : vz,
        "ry" : ry,
        "rz" : rz
    }        
    return velocity