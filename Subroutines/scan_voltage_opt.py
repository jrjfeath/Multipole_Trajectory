import copy
import sys
import time
import numpy as np
import Subroutines.drawing as drawing

from Subroutines.skimmer_trajectory import skimmer_trajectory

def within_boundry(ry,rz,r0):
    '''Calculate if atom/molecule is within specified boundry.'''
    #Radial position
    rint = (ry * ry) + (rz * rz)      
    #Check if molecules within specified boundry
    entered = np.argwhere(rint < (r0 * r0))
    return entered

def within_linear_boundry(time,traj,r=0.001,adjust=False,return_entered=False):
    '''Calculate linear trajectory and see if it is within bounds.'''
    #Calculate off-axis position
    ry = (traj[:,1] * time) + traj[:,3]
    rz = (traj[:,2] * time) + traj[:,4]

    #Check if molecules within specified boundry
    entered = within_boundry(ry,rz,r)

    #Depending on position we may want to overwrite data
    if adjust == True:
        traj[:,3] = ry
        traj[:,4] = rz

    if return_entered == True:
        return traj, entered

    traj = traj[entered][:,0]

    return traj

def within_multipole_boundry(time,traj,sqc,invc,r=0.001,adjust=False,return_entered=False):
    '''Calculate harmonic trajectory and see if it is within bounds.'''
    #Check ry (3) & rz (4) (position) at end of multipole using simple harmonic motion
    ry = (traj[:,3] * np.cos(sqc * time)) + (traj[:,1] * invc * np.sin(sqc * time))
    rz = (traj[:,4] * np.cos(sqc * time)) + (traj[:,2] * invc * np.sin(sqc * time))
    
    #Only keep trajectories that enter
    entered = within_boundry(ry,rz,r)

    #Calculate vy (1) & vz (2) (radial velocity) after exiting multipole
    if adjust == True:
        traj[:,1] = (traj[:,1] * np.cos(sqc * time)) - (traj[:,3] * sqc * np.sin(sqc * time))
        traj[:,2] = (traj[:,2] * np.cos(sqc * time)) - (traj[:,4] * sqc * np.sin(sqc * time))
    
    #Adjust ry (3) & rz (4) (position)
    traj[:,3] = ry
    traj[:,4] = rz

    if return_entered == True:
        return traj, entered

    traj = traj[entered][:,0]

    return traj

def record_data(traj,entered,cpi,gp,positions):
    ry, rz = np.array(traj[:,3]), np.array(traj[:,4])
    gp[entered] = np.add(gp[entered],1)
    track = np.argwhere(gp == cpi+1)[:,0]
    try:
        positions[cpi,track,0] = ry[track]
        positions[cpi,track,1] = rz[track]
    except IndexError:
        print(len(positions),cpi,track)
        sys.exit()
    cpi+=1
    return gp, positions, cpi

def scan_trajectory(d): 
    start = time.time()
    increment = 0.001 # 1mm
    #Generate a 2d and 3d plot classes for plotting our data
    if d['plot_2d_t'] == True: 
        d2 = drawing.d2_drawing()
        d2.find_poles_and_collision(d)
    #Find the distances to draw multipoles and collision area
    if d['plot_3d_t'] == True: 
        d3 = drawing.d3_drawing()
        d3.multipole_radius(d)
    #Create position data (px:(py,pz)) for each trajectory (vest)
    positions = np.empty((int(round(d2.length,3)*1000)+1,d['vest'],2))
    positions[:,:] = np.NAN
    #Good positions
    gp = np.zeros((d['vest'],1))
    #Index of current position
    cpi = 0
    #Index of voltage
    vi = d['V'].index(d['setV'])
    #Calculate trajectories
    traj = skimmer_trajectory(d)

    #Check ry (3) & rz (4) (position) up to skimmer
    r0 = 1 # Set large radius to check against
    increments = np.arange(0,d['lskimmer']+increment,increment)
    for px in increments:
        tohex = px / traj[:,0]
        #When we reach the skimmer check if it passes
        if px == d['lskimmer']: r0 = d['skmr_radius']
        #Overwrite position (cpi) index for all trajectories
        if px != increments[-1]:
            mtraj, entered = within_linear_boundry(tohex,copy.deepcopy(traj),r0,adjust=True,return_entered=True)
            gp, positions, cpi = record_data(mtraj,entered,cpi,gp,positions)
        else:
            traj, entered = within_linear_boundry(tohex,copy.deepcopy(traj),r0,adjust=True,return_entered=True)
            gp, positions, cpi = record_data(traj,entered,cpi,gp,positions)
        #Reset radius after checking collision with skimmer
        if px == d['lskimmer']: r0 = 1

    #Check ry (3) & rz (4) (position) at start of multipole
    increments = np.arange(increment,d['lsource']+increment,increment)
    for px in increments:
        tohex = px / traj[:,0]
        #Overwrite position (cpi) index for all trajectories
        if px != increments[-1]:
            mtraj, entered = within_linear_boundry(tohex,copy.deepcopy(traj),1,adjust=True,return_entered=True)
            gp, positions, cpi = record_data(mtraj,entered,cpi,gp,positions)
        else:
            traj, entered = within_linear_boundry(tohex,copy.deepcopy(traj),1,adjust=True,return_entered=True)
            gp, positions, cpi = record_data(traj,entered,cpi,gp,positions)

    #Loop through each multipole
    for mindex, m in enumerate(d['multipole']):
        #fetch multipole parameters
        p = d['multipole'][m]
        sqc , invc, r0, dist, lpole = p['sqc'], p['invc'],  p['r0'], p['ldist'], p['lpole']
        sqc , invc = sqc[vi], invc[vi]
        
        #Loop through x positions as we travel through multipole
        increments = np.arange(increment,lpole+increment,increment)
        for px in increments:
            thex = px / traj[:,0] # thex = time spent in multipole 
            if px != increments[-1]:
                #If invc == 0 that means the molecule has no dipole and is not effected by the multipole
                if invc == 0:
                    mtraj, entered = within_linear_boundry(thex,copy.deepcopy(traj),r0,adjust=True,return_entered=True)
                else:
                    mtraj, entered = within_multipole_boundry(thex,copy.deepcopy(traj),sqc,invc,r0,adjust=False,return_entered=True)
                gp, positions, cpi = record_data(mtraj,entered,cpi,gp,positions)
            else:
                if invc == 0:
                    traj, entered = within_linear_boundry(thex,copy.deepcopy(traj),r0,adjust=True,return_entered=True)
                else:
                    traj, entered = within_multipole_boundry(thex,copy.deepcopy(traj),sqc,invc,r0,adjust=True,return_entered=True)
                gp, positions, cpi = record_data(traj,entered,cpi,gp,positions)

        #Check if pinhole exists, if true check if it collides
        pin_dist = -1
        if p.get('pin_pos') != None:
            pin_dist = p['pin_pos']

        #Check ry (3) & rz (4) (position) at start of multipole
        r0 = 1
        increments = np.arange(increment,dist+increment,increment)
        for px in increments:
            tohex = px / traj[:,0]
            if px == pin_dist: r0 = p['pin_r']
            #Overwrite position (cpi) index for all trajectories
            if px != increments[-1]:
                mtraj, entered = within_linear_boundry(tohex,copy.deepcopy(traj),r0,adjust=True,return_entered=True)
                gp, positions, cpi = record_data(mtraj,entered,cpi,gp,positions)
            else:
                traj, entered = within_linear_boundry(tohex,copy.deepcopy(traj),r0,adjust=True,return_entered=True)
                gp, positions, cpi = record_data(traj,entered,cpi,gp,positions)
            #After we check if collision into pinhole reset to original multipole radius
            if px == pin_dist: r0 = p['r0']

    #Check ry (3) & rz (4) (position) at end of multipole and to collision region
    increments = np.arange(increment,d['lcollision']+increment,increment)
    for px in increments:
        tohex = px / traj[:,0]
        #Overwrite position (cpi) index for all trajectories
        if px != increments[-1]:
            mtraj, entered = within_linear_boundry(tohex,copy.deepcopy(traj),r0,adjust=True,return_entered=True)
            gp, positions, cpi = record_data(mtraj,entered,cpi,gp,positions)
        else:
            #Collision regions is a 1mm * 1mm area, hence the r0 = 0.001
            traj, entered = within_linear_boundry(tohex,copy.deepcopy(traj),0.0005,adjust=True,return_entered=True)
            gp, positions, cpi = record_data(traj,entered,cpi,gp,positions)

    track = np.argwhere(gp == cpi)[:,0]
    pt = round((len(track) / d["vest"]) * 100,2) # Percent transmission
    print(f'Percent of molecules that make it to collision region: {pt}%')

    print(f'Runtime: {round(time.time() - start,2)}s')

    increments = np.arange(0,d2.length+increment,increment)
    shown = d['traj']
    if d['vest'] < shown: shown = d['vest']
    for i in range(shown):
        rc = np.random.choice((-1,1))
        if d['plot_2d_t'] == True: 
            y = positions[:,i,0]
            z = positions[:,i,1]
            yz = (y ** 2 + z ** 2) ** 0.5
            yz = yz * rc
            d2.draw_molecule(increments, yz, positions[:,i,1])
        if d['plot_3d_t'] == True: 
            d3.draw_molecule(increments, positions[:,i,0], positions[:,i,1])
    y, z = positions[-1,:,0], positions[-1,:,1]
    yz = (y ** 2 + z ** 2) ** 0.5
    yz = yz[np.invert(np.isnan(yz))]
    poo, bins = np.histogram(yz,bins=30)
    bins = bins[:-1] + np.diff(bins)
    hist = drawing.hist_coll()
    hist.draw_hist(bins,poo)

    drawing.draw_plot()

def scan_voltage(d):
    '''
    Scan through a series of voltages to obtain the ideal voltage.\n
    sqc and invc are the derivatives of force applied by the multipole. \n
    r0 is the distance from the centre of the multipole to the rod. \n
    ldist is the distance between the end of the multipole and the next point.
    '''

    start = time.time()

    #Draw an empty plot for data
    plot = drawing.voltage_plot()

    #Calculate an initial trajectory for each of our molecules
    #traj columns are as follows: vel, vy, vz, ry, rz
    traj = skimmer_trajectory(d)

    #Check if trajectory hits skimmer
    toskim = d['lskimmer'] / traj[:,0]
    #Check ry (3) & rz (4) (position) at skimmer
    traj = within_linear_boundry(toskim,traj,d['skmr_radius'],adjust=True)

    #Calculate time taken to reach multipole from skimmer
    tohex = d['lsource'] / traj[:,0]
    #Check ry (3) & rz (4) (position) at start of multipole
    traj = within_linear_boundry(tohex,traj,1,adjust=True)

    entered = len(traj)

    trajs  = [copy.deepcopy(traj) for _ in d['V']]

    #Loop through each multipole
    for mindex, m in enumerate(d['multipole']):
        #fetch multipole parameters
        p = d['multipole'][m]
        sqc , invc, r0, dist, lpole = p['sqc'], p['invc'],  p['r0'], p['ldist'], p['lpole']

        #Calculate position & velocity at entrance of multipole
        #Loop through each set of trajectories belonging to a specific voltage
        for ti, traj in enumerate(trajs):

            #Calculate position & velocity at end of multipole
            thex = lpole / traj[:,0] # thex = time spent in multipole 
            #Check if molecule makes it through multipole
            traj = within_multipole_boundry(thex,traj,sqc[ti], invc[ti],r0,adjust=True)

            #Check if pinhole exists,
            if p.get('pin_pos') != None:
                #Calculate position at distance d from end of multipole
                #We don't overwrite any velocity/positions here as the next multipole checks 
                #from the end of the last multipole; however, we do remove any entries that 
                #hit the area around the pinhole
                tcol = p['pin_pos'] / traj[:,0]
                traj = within_linear_boundry(tcol,traj,p['pin_r'],adjust=False)

            #Distance to a subsequent multipole
            if dist != 0:
                tcol = dist / traj[:,0]
                traj = within_linear_boundry(tcol,traj,1,adjust=True)

            trajs[ti] = traj

    for ti, traj in enumerate(trajs):
        tcol = d['lcollision'] / traj[:,0]
        traj = within_linear_boundry(tcol,traj)
        trajs[ti] = traj

    print(time.time()-start)

    #Plot the data
    x = [i for i in d['V']]
    y = [(len(i) / entered) * 100 for i in trajs]
    plot.draw(x,y)
    for i, v in enumerate(x):
        print(v,y[i])

    drawing.draw_plot()  
    #Destroy the plot in the event user runs multiple scan commands
    plot.destroy()
    