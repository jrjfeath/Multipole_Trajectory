import copy
import time
import numpy as np
import Subroutines.drawing as drawing

from Subroutines.skimmer_trajectory import skimmer_trajectory

#Speed of light
c = 299792458  #m * s_1
#Plank's constant
h = 6.626068e-34 #m^2 * kg * s_1
hbar = 1.05457148e-34 #m^2 * kg * s_1 = (h/2pi)
hc = h * c #m^3 * kg * s_2

def outside_boundry(rp,r0):
    '''Calculate if atom/molecule is outside specified boundry.'''    
    #Check if molecules outside specified boundry
    return np.argwhere(abs(rp) >= r0)

def within_boundry(rp,r0):
    '''Calculate if atom/molecule is within specified boundry.'''    
    #Check if molecules outside specified boundry
    return np.argwhere(abs(rp) <= r0)

def within_linear_boundry(time,traj,r=0.1):
    '''Calculate linear trajectory and see if it is within bounds.'''
    #Calculate radial position
    traj[:,2] = (traj[:,1] * time) + traj[:,2]
    ob = outside_boundry(traj[:,2], r)
    traj[ob] = np.array([np.nan,np.nan,np.nan])
    return traj

def within_multipole_boundry(t,t2,vf1,vf2,traj,mass,r=0.001):
    '''Calculate radial trajectory in multipole and see if it is within bounds.'''
    #Calculate the radial acceleration of the molecule
    ar = ((traj[:,2] * vf1) * (1 / (1 + (vf2 / traj[:,2]**2)**2)**0.5)) / mass
    #Calculate the new radial position of the molecule
    traj[:,2] += ((ar * t2) / 2) + (traj[:,1] * t)
    #Calculate the new radial velocity of the molecule
    traj[:,1] += (ar * t)
    #Calculate which trajectories are outside the multipole
    ob = outside_boundry(traj[:,2], r)
    traj[ob] = np.array([np.nan,np.nan,np.nan])
    return traj

def scan_voltage(d):
    '''
    Scan through a series of voltages to obtain the ideal voltage.\n
    sqc and invc are the derivatives of force applied by the multipole. \n
    r0 is the distance from the centre of the multipole to the rod. \n
    ldist is the distance between the end of the multipole and the next point.
    '''
    #Generate a 2d and 3d plot classes for plotting our data
    if d['plot_2d_t'] == True: 
        d2 = drawing.d2_drawing()
        d2.find_poles_and_collision(d)
    #Find the distances to draw multipoles and collision area
    if d['plot_3d_t'] == True: 
        d3 = drawing.d3_drawing()
        d3.multipole_radius(d)

    start = time.time()
    increment = 0.010 # 5mm

    #Calculate an initial trajectory for each of our molecules
    #traj columns are as follows: longitudinal velocity, radial velocity, radial position
    traj = skimmer_trajectory(d)

    x = [[0 for _ in range(len(traj[:,2]))]]
    y = [np.copy(traj[:,2])]

    #Check if trajectory hits skimmer
    toskim = d['lskimmer'] / traj[:,0]
    #Check radial position at skimmer
    traj = within_linear_boundry(toskim,traj,r=d['skmr_radius'])

    x.append([d['lskimmer'] for _ in range(len(traj[:,2]))])
    y.append(np.copy(traj[:,2]))

    #Calculate time taken to reach multipole from skimmer
    tohex = d['lsource'] / traj[:,0]
    #Check radial position at start of multipole
    traj = within_linear_boundry(tohex,traj)

    start = d['lskimmer']+d['lsource']

    x.append([start for _ in range(len(traj[:,2]))])
    y.append(np.copy(traj[:,2]))

    #Make a copy of trajectories for each voltage 
    trajs  = [copy.deepcopy(traj) for _ in d['V']]

    #Loop through each multipole and determine the effect on the trajectory
    for key in d['multipole']:
        m = d['multipole'][key]
        #Loop through voltages to determine effect on molecule
        for vi in range(len(d['V'])):
            traj = trajs[vi]
            vf1 = m['f1'][vi]
            vf2 = m['f2'][vi]
            #Calculate time spent in hexapole so far
            t = increment / traj[:,0]
            t2 = t * t
            #Loop through x positions as we travel through multipole
            increments = np.arange(increment,m['lpole']+increment,increment)
            for px in increments:
                traj = within_multipole_boundry(t, t2, vf1, vf2, traj, d['mass'], m['r0'])
                #If the user is tracking trajectory of a voltage
                if (vi == d['Vindex']) or (len(m['f1']) == 1):
                    x.append([start+px for _ in range(len(traj[:,2]))])
                    y.append(np.copy(traj[:,2]))

            #Check if pinhole exists
            if m.get('pin_pos'):
                tcol = m['pin_pos'] / traj[:,0]
                traj = within_linear_boundry(tcol,traj)
                #If the user is tracking trajectory of a voltage
                if (vi == d['Vindex']) or (len(m['f1']) == 1):
                    x.append([start+m['lpole']+m['pin_pos'] for _ in range(len(traj[:,2]))])
                    y.append(np.copy(traj[:,2]))
                #Determine which trajectories hit collision region
                if m['pin_r'] > 0:
                    ob = outside_boundry(traj[:,2],m['pin_r'])
                else:
                    ob = within_boundry(traj[:,2],abs(m['pin_r']))
                traj[ob] = np.array([np.nan,np.nan,np.nan])

            #Check if there is a second multipole and the distance to it
            if m.get('ldist'):
                dist = m['ldist']
                #If there is a pinhole subtract the distance calculated for that
                if m.get('pin_pos'):
                    tcol = (dist-m['pin_pos']) / traj[:,0]
                else:
                    tcol = dist / traj[:,0]
                traj = within_linear_boundry(tcol,traj)
                #If the user is tracking trajectory of a voltage
                if (vi == d['Vindex']) or (len(m['f1']) == 1):
                    x.append([start+m['lpole']+dist for _ in range(len(traj[:,2]))])
                    y.append(np.copy(traj[:,2]))     

            #If vset exists set all trajectories to match this voltage
            if m.get('vset'):
                trajs  = [copy.deepcopy(traj) for _ in d['V']]
                break

            #Record the trajectory at the end of the multipole & pinhole
            trajs[vi] = traj

        start+=m['lpole']
        if m.get('ldist'): start+=dist

    start+=d['lcollision']
    counts = []
    for ti, traj in enumerate(trajs):
        tcol = d['lcollision'] / traj[:,0]
        traj = within_linear_boundry(tcol,traj)
        #Draw before we remove trajectories that don't hit collision regions
        if ti == d['Vindex']:
            x.append([start for _ in range(len(traj[:,2]))])
            y.append(np.copy(traj[:,2]))
            #Plot histogram of rp at collision region
            bins = np.arange(-d['crr'],d['crr']+d['crr']/10,d['crr']/10)
            d2.draw_hist(bins,traj[:,2])
        #Determine which trajectories hit collision region
        ob = outside_boundry(traj[:,2],d['crr'])
        traj[ob] = np.array([np.nan,np.nan,np.nan])
        transmission = (d['vest']-sum(np.isnan(traj[:,2])))/d['vest']
        counts.append(transmission * 100)

    if d['calcv'] == True:
        #Draw an empty plot for data
        plot = drawing.voltage_plot()
        #Plot transmission data
        plot.draw(d['V'],counts)
    
    if d['plot_2d_t'] == True:
        #Draw trajectories on plot
        d2.draw_molecule(np.array(x),np.array(y))

    if d['plot_3d_t'] == True:
        pass
        # d3.draw_molecule(x,y,y)
    
    drawing.draw_plot()