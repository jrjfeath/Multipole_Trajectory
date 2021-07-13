import copy
import multiprocessing
import os
import random
import sys
import time
import yaml

from math import sqrt, cos, sin, pi, log, exp

import matplotlib.pyplot as plt
import numpy as np
import psutil

def calculate_c(d,n,r0):
    '''Calculate the force applied to the molecule done by the multipole.\n
       See additional info for more details into the method.
    '''
    #Pull variables from dict, faster than calling every loop
    dip, mass, j = d['dipole'], d['mass'], d['j']
    m, epi, _lambda = d['m'], d['epi'], d['lambda']
    #Make empty lists for saving results
    sqc, invc = [], []
    for V in d['V']:
        c = (dip * n * (n-1) * V) / (mass * r0**n) * (abs(epi * m *_lambda) / (j * (j+1)))
        #Calculate some derivatives of c to speed up calculations
        sqc.append(sqrt(c))
        invc.append(1/sqc[-1])
    return sqc, invc

def positive_or_negative():
    '''Returns a random positive or negative.'''
    return 1 if random.random() < 0.5 else -1

def skimmer_trajectory(d):
    '''Determine the trajectory of molecules out of the skimmer.'''
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
    #Calculate the off-axis positions if requested
    if d['srcd'].strip().lower() == 's': #single distribution
        ry = d['fwhmsrc']
        rz = d['fwhmsrc']
    elif d['srcd'].strip().lower() == 'u': #uniform distribution
        r = positive_or_negative() #Grab a random +/-
        ry = d['fwhmsrc'] * urnd * r
        rz = d['fwhmsrc'] * urnd * r 
    else: #gaussian distribution
        ry = d['fwhmsrc']*grnd
        rz = d['fwhmsrc']*grnd           
    return vel,vy,vz,ry,rz  

def check_radial_position(r0,ry,rz):
    '''Check if molecule is within bounds.\n
       r0 = radius of the multipole\n
       ry,rz = off-axis position of molecule
    '''
    rint = (rz * rz) + (ry * ry)
    #if molecule is outside bounds return true
    if(rint >= (r0 * r0)): return 1
    return 0

def multipole_trajectory(dist,lpole,vel,vy,vz,ry,rz,sqc,invc,r0):
    '''Determine if molecule enters and exits the multipole successfully.\n
       dist = distance from skimmer or multipole to entrance of multipole\n
       lpole = length of the multipole\n
       vel = velocity of molecule\n
       vy,vz = off-axis velocities\n
       ryp,ryz = off-axis positions\n
       sqc,invc = functions of c precalculated to save time
       r0 = radius of the multipole in meters
    '''
    info = {'enter':0,'exit':0,'ryhex':0,'rzhex':0,'velyf':0,'velzf':0}
    #Calculate position & velocity at entrance of hexapole
    tohex = dist/vel #Calculate time taken to reach hexapole
    ry0 = (vy * tohex) + ry
    rz0 = (vz * tohex) + rz
    #Check if the molecule will enter the hexapole
    if check_radial_position(r0,ry0,rz0) == 1: return info
    info['enter']=1
    
    #Calculate position & velocity at end of hexapole
    thex = lpole/vel # thex = time spent in hexapole 
    ryhex = (ry0 * cos(sqc * thex)) + (vy * invc * sin(sqc * thex))
    rzhex = (rz0 * cos(sqc * thex)) + (vz * invc * sin(sqc * thex))
    #Check if the molecule will exit the hexapole
    if check_radial_position(r0,ryhex,rzhex) == 1: return info
    info['exit']=1

    velzf = (vz * cos(sqc * thex)) - (rz0 * sqc * sin(sqc * thex))
    velyf = (vy * cos(sqc * thex)) - (ry0 * sqc * sin(sqc * thex))
    
    info['ryhex'] = ryhex
    info['rzhex'] = rzhex
    info['velyf'] = velyf
    info['velzf'] = velzf

    return info

def collision_trajectory(vel,vy,vz,ry,rz,lcol):
    '''Calculate if molecule hits collision area.'''
    #Calculate position at distance d from end of hexapole & increment image intensity
    tcol = lcol / vel
    ryf = (vy * tcol) + ry
    rzf = (vz * tcol) + rz   

    # is the trajectory in the interaction region ( 1mm x 2mm x 1mm ) --> 2mm in beam velocity direction
    if ( ( rzf <= 0.001 ) and (rzf >= -0.001 ) and ( ryf <= 0.001 ) and ( ryf >= -0.001 ) ): 
        return 1
    return 0

def search_parameters(d,densities):
    '''Calculate the number of molecules that reach each point.'''
    #Calculate an initial trajectory for each of our molecules
    for _ in range(d['est']):
        #What did densities look like before this trajectory?
        initial = copy.deepcopy(densities)
        lsource = d['lsource']
        vel, vy, vz, ry, rz = skimmer_trajectory(d)
        vra = np.zeros((len(d['V']),4))
        #Check if molecule makes it through each multipole
        for column, m in enumerate(d['multipole']): 
            p = d['multipole'][m] #multipole parameters
            sqc , invc = d['multipole'][m]['sqc'], d['multipole'][m]['invc']
            r0 = d['multipole'][m]['r0']
            #Change the voltage on the multipoles
            for row in range(len(d['V'])):
                cid = (column * 2)-1
                #Check if we're on the first multipole
                if column == 0:
                    info = multipole_trajectory(lsource,p['lpole'],vel,vy,vz,ry,rz,sqc[row],invc[row],r0)
                #Check if molecule did not successfully exit previous multipole
                elif column != 0 and (densities[row][cid] - initial[row][cid]) == 0:
                    continue
                else: #If it successfully passed, calculate trajectory
                    #We need to use the velocities and positions from the previous multipole and
                    #not the source
                    vmy, vmz, rmy, rmz = list(vra[row])
                    info = multipole_trajectory(lsource,p['lpole'],vel,vmy,vmz,rmy,rmz,sqc[row],invc[row],r0)
                #overwrite velocity and position data
                vra[row] = [info['velyf'],info['velzf'],info['ryhex'],info['rzhex']]
                #write success to densities
                densities[row][column] += info['enter']
                densities[row][column+1] += info['exit']
            #Change the source distance
            lsource = d['multipole'][m]['ldist']
        #Calculate if molecule hits collision area
        for row in range(len(d['V'])):
            if densities[row][-2] - initial[row][-2] == 0: continue
            vmy, vmz, rmy, rmz = list(vra[row])
            hit = collision_trajectory(vel,vmy,vmz,rmy,rmz,d['lcollision'])
            densities[row][-1]+=hit
    collides = max(densities[::,2])
    index = np.where(densities[:,2] == collides)
    row = densities[index]
    V = [d['V'][x] for x in index[0]]
    return densities,row,V

def plot_data(x,y,xl,yl,title):
    _, ax = plt.subplots()
    ax.plot(x,y)
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ax.set_title(title)
    plt.show()

def load_input():
    '''Load parameter file containing user input.'''
    #Get directory containing files
    cfd = os.path.dirname(__file__)
    #Before we proceed check if we can load Parameter file
    try:
        with open(os.path.join(cfd,'Parameters.yaml'), 'r') as stream:
            d = yaml.load(stream,Loader=yaml.SafeLoader)
    except FileNotFoundError:
        print("Cannot find the Parameters file, where did it go?")
        sys.exit()
    #If file exists modify input values
    # convert to SI unit 
    atomic_mass = 1.6605390666e-27 #au
    debye = 3.33654e-30 #Cm, electric dipole moment
    d['mass']*=atomic_mass
    d['dipole']*=debye
    d['fwhmsk']*=(pi / 180.0) #convert to radians
    d['tpi'] = pi*2 #calculate pi * 2 to save time
    #Convert lengths into meters
    d['lsource']/=1000
    d['lcollision']/=1000
    d['lmax']/=1000
    d['sd']/=1000
    return d

def setup():
    '''Load user parameters and modify as needed.'''
    #Load the parameters specified by the user
    d = load_input()
    #Check if user wants to scan voltages
    if d['calcv'] == True: 
        d['V'] = range(d['Vstart'],d['Vfinal']+1,d['Vstep'])
    else: 
        d['V'] = [d['setV']]
    counts = 1 #How many densities to keep track of?
    #Precalculate the values of the force applied by multipole (c)
    for key in d['multipole']:
        m = d['multipole'][key]
        n = m['size'] / 2 #Half the size of the multipole, i.e n = 3 for hexapole
        r0 = (m['d0'] / 2) / 1000 #Radius of the multipole
        d['multipole'][key]['r0'] = r0
        sqc, invc = calculate_c(d,n,r0)
        d['multipole'][key]['sqc'] = sqc
        d['multipole'][key]['invc'] = invc
        #Convert lengths to meters
        d['multipole'][key]['lpole']/=1000
        d['multipole'][key]['ldist']/=1000
        counts+=2
    d['counts'] = counts
    return d

def randomize_distances(d):
    md = d['lmax'] #max distance
    for key in d['multipole']:
        d['multipole'][key]['lpole']
    return d

if __name__ == "__main__":
    d = setup()
    #Make a list of the number of molecules at each position
    densities = np.zeros((len(d['V']),d['counts']))
    if d['search'] == False:
        start = time.time()
        densities,row,V = search_parameters(d,copy.deepcopy(densities))
        finish = round(time.time()-start,2)
        #Output best voltages, transmission data for voltages, and runtime
        print(V,row,finish)
        x = [x for x in d['V']]
        y = (densities[::,2] / d['est']) * 100
        xl = 'Voltages (V)'
        yl = 'Transmission (%)'
        title = 'What percent of molecules are transmitted to collision region?'
        plot_data(x,y,xl,yl,title)
    else:
        e = randomize_distances(copy.deepcopy(d))
        '''#Setup multiprocessing to speed up process
        q_in = multiprocessing.Queue()
        q_out = multiprocessing.Queue()

        #Only allow as many workers as there are cores
        nproc = d['nproc']
        if nprocs > psutil.cpu_count():
            nprocs = psutil.cpu_count()        
        
        proc = [multiprocessing.Process(target=Opt_Input, args=(hdf_data, q_in, q_out)) for i in range(nprocs)]'''