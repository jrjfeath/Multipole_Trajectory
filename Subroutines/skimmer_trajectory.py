import numpy as np

def skimmer_trajectory(d):
    '''
    Determine the trajectory of molecules out of the skimmer.\n
    Returns a dictionary containing the velocity in the x, y, and z directions.\n
    Additionally returns the off axis position of the molecules.\n
    traj = [vx,vy,vz,ry,rz]
    '''
    #Calculate a random gaussian velocity based around theoretical velocity
    grnd = np.random.normal(0, 1, d['vest']) #gaussian random +/- 2.5
    urnd = np.random.uniform(0, 1, d['vest']) #uniform random distribution
    grnda = grnd * d['fwhmsk'] #Calculate angle of molecules
    vel = d['velocity'] + (grnd * d['fwhmv']) #Calculate velocity towards collision
    
    #Calculate a uniform float between 0,1 to simulate an even distribution
    #of velocity angles, multiply by 2pi to get -/+ angles
    urnda = urnd * d['tpi']
    #Calculate velocities in the y & z directions
    vy = vel * np.sin(grnda) * np.sin(urnda) 
    vz = vel * np.sin(grnda) * np.cos(urnda)
    
    cy, cz = d.get('skmr_pos') #Grab off-axis positions
    dist = d.get('skmr_dist') #Check what kind of distribution was selected
    skmr_r = d.get('skmr_radius') * 1000 #Grab skimmer radius

    #If no distribution is select centre around specified position
    if dist == 's':
        ry = np.zeros((1,d['vest']))[0]
        ry.fill(cy)
        rz = np.zeros((1,d['vest']))[0]
        rz.fill(cz)

    #Calculate a Gaussian centred at 0 with 90% of the integral in [-sr,sr]
    #We need to divide the radius by 2 to get even proportions in -/+
    if dist == 'g':
        sigma = (skmr_r / 2 / 1.644854)
        ry = np.random.normal(cy, sigma, d['vest']) / 1000
        rz = np.random.normal(cz, sigma, d['vest']) / 1000

    #Calculate a uniform distribution with limits of skmr_r
    if dist == 'u':
        ry = cy + np.random.uniform(-skmr_r, skmr_r, d['vest']) / 1000
        rz = cz + np.random.uniform(-skmr_r, skmr_r, d['vest']) / 1000

    traj = np.stack((vel,vy,vz,ry,rz), axis=-1)

    return traj