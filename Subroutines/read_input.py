import math
import os
import sys
import yaml
import Subroutines.forces as forces
import numpy as np

def setup(filename):
    '''
    Load parameter file containing user input.
    '''
    #Before we proceed check if we can load Parameter file
    try:
        with open(filename, 'r') as stream:
            d = yaml.load(stream,Loader=yaml.SafeLoader)
    except FileNotFoundError:
        print("Cannot find the Parameters file, where did it go?")
        sys.exit()

    d['filename'] = filename
    d['Stark_Calculated'] = False #Have the stark effects been written to file?
    
    # variables to convert to SI unit 
    atomic_mass = 1.6605390666e-27 #(kg)
    debye = 3.33654e-30 #Cm, electric dipole moment
    ec = 1.602176462e-19 # electron charge (C/e)
    #Speed of light
    c = 299792458  #m * s_1
    #Plank's constant
    h = 6.626068e-34 #m^2 * kg * s_1
    hc = c * h
    d['tpi'] = 2 * math.pi #frequency per cycle(rad/cycle)
    
    #Apply SI conversions
    d['mass']*=atomic_mass
    d['dipole']*=debye
    d['charge']*=ec
    d['Vfrequency'] = float(d['Vfrequency']) * d['tpi']

    #convert to radians
    d['fwhmsk']*=(math.pi / 180.0) 

    #Convert lengths into meters
    d['skmr_radius']/=1000
    d['lskimmer']/=1000
    d['lsource']/=1000
    d['lcollision']/=1000

    #Check if molecule is linear
    if d['k'] == 0: 
        epi, omega = d['epi'], d['omega']
        d['k'] = epi * omega

    #Calculate kappa value
    d['kappa'] = (d['m'] * d['k']) / (d['j'] * (d['j'] + 1))
    #Check if the user has specified lambda doublet splitting
    if 'w_inv' not in d: d['w_inv'] = 0
    d['w_inv'] *= (hc * 100) # In Joules

    #Check if user wants to scan voltages, if true make a list containing voltages
    if d['calcv'] == True: 
        d['V'] = np.array(range(d['Vstart'],d['Vfinal']+1,d['Vstep']))
    else: 
        d['V'] = [d['setV']]

    #Determine and modify multipole variables
    counts = 1 #How many densities to keep track of?
    #Precalculate the values of the force applied by multipole (c)
    #key corresponds to each respective multipole
    dip, mass = d['dipole'], d['mass']
    for key in d['multipole']:
        #fetch multipole parameters
        m = d['multipole'][key]
        #Convert diameter to m and then divide by 2 for radius
        r0, size = m['d0'] / (2 * 1000), m['size']
        #Calculate the unchanging parts of the acceleration equation for multipoles
        m['f1'] = (size * dip * d['kappa'] * d['V']) / (mass * r0 * r0 * r0)
        f2 = (d['w_inv'] * r0 * r0 * r0) / (size * dip * d['kappa'] * d['V'])
        #If the user checks 0 volts then f2 will equal negative infinity
        f2[f2 == -np.inf] = 0
        m['f2'] = f2

        #Convert lengths to meters
        m['lpole']/=1000
        m['ldist']/=1000
        if 'pin_pos' in m: m['pin_pos']/=1000
        if 'pin_r' in m: m['pin_r']/=1000
        if 'pinhole' in m: counts+=1
        counts+=2
    d['counts'] = counts

    return d