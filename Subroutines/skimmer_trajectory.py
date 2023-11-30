import numpy as np

def skimmer_trajectory(d):
    #Calculate a random gaussian velocity based around theoretical velocity
    grnd = np.random.normal(0, 1, d['vest']) #gaussian random +/- 2.5
    grnda = grnd * d['fwhmsk'] #Calculate angle of molecules
    vel = d['velocity'] + (grnd * d['fwhmv']) #Calculate total velocity
    rv = vel * np.sin(grnda) #Calculate radial velocity
    lv = vel * np.cos(grnda) #Calculated longitudinal velocity

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
    #Calculate the radial position
    rp = (ry ** 2 + rz ** 2)**0.5

    traj = np.stack((lv,rv,rp), axis=-1)

    return traj