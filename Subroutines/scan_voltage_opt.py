import numpy as np
import Subroutines.drawing as drawing

from Subroutines.check_collision import check_collision
from Subroutines.multipole_trajectory import multipole_trajectory
from Subroutines.skimmer_trajectory import skimmer_trajectory

def scan_voltage(d):
    '''
    Scan through a series of voltages to obtain the ideal voltage.\n
    sqc and invc are the derivatives of force applied by the multipole. \n
    r0 is the distance from the centre of the multipole to the rod. \n
    ldist is the distance between the end of the multipole and the next point.
    '''

    #Draw an empty plot for data
    plot = drawing.voltage_plot()

    #Make an array to hold the data at each voltage and each step
    densities = np.zeros((len(d['V']),d['counts']))

    #Calculate an initial trajectory for each of our molecules
    #traj columns are as follows: vel, vy, vz, ry, rz
    traj = skimmer_trajectory(d)

    #Loop through each multipole
    for mindex, m in enumerate(d['multipole']):
        #fetch multipole parameters
        p = d['multipole'][m]
        sqc , invc, r0, dist, lpole = p['sqc'], p['invc'],  p['r0'], p['ldist'], p['lpole']
        
        #fetch the initial distance from skimmer to multipole
        if mindex == 0:
            dist = d['lsource']

        #Calculate position & velocity at entrance of multipole
        tohex = dist/traj[:,0] #Calculate time taken to reach multipole
        traj[:,3] = (traj[:,1] * tohex) + traj[:,3]
        traj[:,4] = (traj[:,2] * tohex) + traj[:,4]
        rint = (traj[:,3] * traj[:,3]) + (traj[:,4] * traj[:,4])
        #Check if molecules enter the multipole
        entered = np.argwhere(rint < (r0 * r0))
        #Only keep trajectories that enter
        traj = traj[entered][:,0]

        #Calculate position & velocity at end of hexapole
        thex = lpole/traj[:,0] # thex = time spent in hexapole 
        #Calculate the trajectory of a molecule/atom through the multipole
        #using simple harmonic motion
        sqc = np.reshape(sqc,(len(sqc),1))
        invc = np.reshape(invc,(len(invc),1))
        ryhex = (traj[:,3] * np.cos(sqc * thex)) + (traj[:,1] * invc * np.sin(sqc * thex))
        rzhex = (traj[:,4] * np.cos(sqc * thex)) + (traj[:,2] * invc * np.sin(sqc * thex))
        #Check if the molecule will exit the hexapole
        print(ryhex)