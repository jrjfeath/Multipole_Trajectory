import Subroutines.check_radial as check_radial
from math import sin, cos

def multipole_trajectory(dist,lpole,vel,vy,vz,ry,rz,sqc,invc,r0):
    '''
    Determine if molecule enters and exits the multipole successfully.\n
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
    if check_radial.check_radial_position(r0,ry0,rz0) == 1: return info
    info['enter']=1
    
    #Calculate position & velocity at end of hexapole
    thex = lpole/vel # thex = time spent in hexapole 
    #Calculate the trajectory of a molecule/atom through the multipole
    #using simple harmonic motion
    ryhex = (ry0 * cos(sqc * thex)) + (vy * invc * sin(sqc * thex))
    rzhex = (rz0 * cos(sqc * thex)) + (vz * invc * sin(sqc * thex))
    #Check if the molecule will exit the hexapole
    if check_radial.check_radial_position(r0,ryhex,rzhex) == 1: return info
    info['exit']=1

    velzf = (vz * cos(sqc * thex)) - (rz0 * sqc * sin(sqc * thex))
    velyf = (vy * cos(sqc * thex)) - (ry0 * sqc * sin(sqc * thex))
    
    info['ryhex'] = ryhex
    info['rzhex'] = rzhex
    info['velyf'] = velyf
    info['velzf'] = velzf

    return info