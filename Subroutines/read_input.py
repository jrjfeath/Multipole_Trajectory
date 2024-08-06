import math
import os
import sys
import yaml
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

    #convert to radians
    d['skmr_std']*=(math.pi / 180.0) 

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
    d['kappa'] = d['m'] * d['k'] / (d['j'] * (d['j'] + 1))
    #Check if the user has specified lambda doublet splitting
    if 'w_inv' not in d: d['w_inv'] = 0
    d['w_inv'] *= (hc * 100) # In Joules

    #Create lists for scanning voltage
    if d.get('setV'): 
        d['V'] = [d['setV']]
    else:
        d['V'] = [0]
    d['Vindex'] = 0
    #If the user is scanning a set of voltages
    if d['calcv'] == True: 
        d['V'] = np.array(range(d['Vstart'],d['Vfinal']+1,d['Vstep']))
    #If the user is checking the trajectories of a specific voltage
    if d['calct'] == True:
        try: d['Vindex'] = np.where(d['V'] == d['setV'])[0][0]
        except IndexError: 
            d['V'] = np.sort(np.append(d['V'],d['setV']))
            d['Vindex'] = np.where(d['V'] == d['setV'])[0][0]

    #Determine and modify multipole variables
    counts = 1 #How many densities to keep track of?
    #key corresponds to each respective multipole
    dip, mass = d['dipole'], d['mass']
    for key in d['multipole']:
        #fetch multipole parameters
        m = d['multipole'][key]
        #Convert diameter to m and then divide by 2 for radius
        r0, size = m['d0'] / (2 * 1000), m['size']
        d['multipole'][key]['r0'] = r0
        if m.get('vset'): v_r = np.array([m['vset']])
        else: v_r = d['V']
        #Calculate the unchanging parts of the acceleration equation for multipoles
        m['f1'] = ((size * dip * d['epi'] * abs(d['kappa']) * v_r) / (r0 * r0 * r0)) 
        f2 = (d['w_inv'] * r0 * r0 * r0) / (size * dip * d['kappa'] * v_r)
        #If the user checks 0 volts then f2 will equal negative infinity
        f2[f2 == -np.inf] = 0
        m['f2'] = f2 

        #Convert lengths to meters
        m['lpole']/=1000
        if 'ldist' in m: m['ldist']/=1000
        if 'pin_pos' in m: m['pin_pos']/=1000
        if 'pin_r' in m: m['pin_r']/=1000
        if 'pinhole' in m: counts+=1
        counts+=2
    d['counts'] = counts

    return d