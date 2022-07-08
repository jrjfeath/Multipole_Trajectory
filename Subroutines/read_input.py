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

    #Check if user wants to scan voltages, if true make a list containing voltages
    if d['calcv'] == True: 
        d['V'] = range(d['Vstart'],d['Vfinal']+1,d['Vstep'])
    else: 
        d['V'] = [d['setV']]

    #Determine and modify multipole variables
    counts = 1 #How many densities to keep track of?
    #Precalculate the values of the force applied by multipole (c)
    #key corresponds to each respective multipole
    for key in d['multipole']:
        m = d['multipole'][key]
        n = m['size'] / 2 #Half the size of the multipole, i.e n = 3 for hexapole
        r0 = (m['d0'] / 2) / 1000 #Radius of the multipole
        m['r0'] = r0

        #Determine the force of the multipole on the species, check if molecule or ion
        if d['charge'] != 0: 
            sqc, invc = forces.calculate_d(d,n,r0)
        else: 
            #Calculate the first and second stark effects and add them
            w_e1 = forces.first_order(d,n,r0)
            w_e2 = forces.second_order(d,n,r0)
            
            #Add the stark effects together
            w = {}
            sqc, invc = [], []
            for E in d['V']:
                w[E] = w_e1[E] + w_e2[E]
                #To apply multipole effects introduce mass, radius, and size (n) = s_2
                #Get units of s_1 (in this case E = U so units of m^3 * kg * s_2)
                #See Ref. 1 for additional details (pg. 7666 - 7667)
                s_1 = np.sqrt(w[E] * hc * ((n * (n-1)) / (d['mass'] * r0 * r0 * r0)))
                sqc.append(s_1)
                #get units of s by taking the inverse of s_1
                if s_1 == 0.0:
                    s = 0
                else:
                    s = 1 / s_1
                invc.append(s)
            
            if d['Stark_Calculated'] == False:
                #Write to the string holding stark data
                strg = 'Voltage, 1st Stark (cm-1), 2nd Stark (cm-1)\n'
                for key in w_e1.keys():
                    strg += '{:05d},{:06f},{:06f}\n'.format(key,w_e1[key],w_e2[key])
                
                #Write stark data to stark.csv
                with open(f'{d["filename"][:-5]}-Stark.csv','w') as opf:
                    opf.write(strg)
                d['Stark_Calculated'] = True

        m['sqc'] = sqc
        m['invc'] = invc

        #Convert lengths to meters
        m['lpole']/=1000
        m['ldist']/=1000
        if m.get('pin_pos') != None: m['pin_pos']/=1000
        if m.get('pin_r') != None: m['pin_r']/=1000
        if m.get('pinhole') == True: counts+=1
        counts+=2
    d['counts'] = counts

    return d