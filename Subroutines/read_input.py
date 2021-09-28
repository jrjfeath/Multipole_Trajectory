import math
import os
import sys
import yaml
import Subroutines.forces as forces

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
            w_e1, sqc1, invc1 = forces.first_order(d,n,r0)
            w_e2, sqc2, invc2 = forces.second_order(d,n,r0)
            
            #Add the effects together
            sqc = sqc1 + sqc2
            invc = invc1 + invc2
            
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