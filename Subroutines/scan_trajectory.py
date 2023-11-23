import numpy as np
import Subroutines.drawing as drawing

from math import cos, sin
from Subroutines.check_radial import check_radial_position
from Subroutines.skimmer_trajectory import skimmer_trajectory

def linear_trajectory(velocity,distance,xd,yd,zd):
    #Check if after skimmer
    try: 
        d0 = xd[-1]
        ry0 = yd[-1]
        rz0 = zd[-1]
    #If before skimmer set defaults
    except IndexError: 
        d0 = 0
        ry0 = velocity['ry']
        rz0 = velocity['rz']
    for d in np.arange(d0,distance,0.001):
        tohex = d / velocity['vx'] 
        ry = (velocity['vy'] * tohex) + ry0
        rz = (velocity['vz'] * tohex) + rz0
        xd.append(d+d0)
        yd.append(ry)
        zd.append(rz)
    return xd, yd, zd

def scan_trajectory(md):

    #Generate a 2d and 3d plot classes for plotting our data
    d2 = drawing.d2_drawing()
    d3 = drawing.d3_drawing()

    #Find the distances to draw multipoles and collision area
    d2.find_poles_and_collision(md)
    #Find the maximum radius of the multipoles for 3d drawing
    d3.multipole_radius(md)

    #get the voltage index to be examined
    vi = md['V'].index(md['setV'])

    est = md['traj']
    count = 0
    skimmer = 0
    for i in range(est):
        #Make empty lists to store our x,y,z data
        xd, yd, zd = [], [], []

        #Calculate an initial trajectory for each of our molecule
        velocity = skimmer_trajectory(md)

        xd, yd, zd = linear_trajectory(velocity,md['lskimmer']+0.001,xd,yd,zd)

        if check_radial_position(md['skmr_radius'],yd[-1],zd[-1]) == 1: 
            d2.draw_molecule(xd, yd, zd)
            d3.draw_molecule(xd, yd, zd)
            skimmer+=1
            continue

        #Calculate position just before the multipole using velocity 
        #out of the skimmer
        xd, yd, zd = linear_trajectory(velocity,md['lsource']+0.001,xd,yd,zd)
        
        #Check if the molecule crashes at some point in time
        crashed = False

        #Loop through our multipoles
        for mindex, m in enumerate(md['multipole']):
            if crashed == True: break

            p = md['multipole'][m] #fetch the multipole parameters

            #Check if molecule enters the multipole
            if check_radial_position(p['r0'],yd[-1],zd[-1]) == 1: break

            d0 = xd[-1] #where does the multipole start?           

            #derivatives of forces applied at voltage V (setV)
            sqc, invc = p['sqc'][vi], p['invc'][vi]
            
            #grab entry multipole positions in the y and z directions
            ry0, rz0 = yd[-1], zd[-1]
            for d in np.arange(0,p['lpole']+0.001,0.001):
                # thex = time spent in multipole
                thex = d / velocity['vx'] 

                ry = (ry0 * cos(sqc * thex)) + (velocity['vy'] * invc * sin(sqc * thex))
                rz = (rz0 * cos(sqc * thex)) + (velocity['vz'] * invc * sin(sqc * thex))
                
                xd.append(d0+d), yd.append(ry), zd.append(rz)

                #Check if molecule hits the multpole
                if check_radial_position(p['r0'],yd[-1],zd[-1]) == 1: 
                    crashed = True
                    break

            #If the molecule crashed in the multipole exit loop
            if crashed == True: break

            #Update d0 to distance after multipole
            d0 = xd[-1]
            
            #Calculate the molecule velocities leaving the multipole
            velzf = (velocity['vz'] * cos(sqc * thex)) - (rz0 * sqc * sin(sqc * thex))
            velyf = (velocity['vy'] * cos(sqc * thex)) - (ry0 * sqc * sin(sqc * thex))

            #Check if we have another loop to do, if not calculate collision
            if mindex == len(md['multipole']) - 1: length = md['lcollision']
            else: length = p['ldist']
            
            #grab exit multipole positions in the y and z directions
            ry0, rz0 = yd[-1], zd[-1]
            for d in np.arange(0,length+0.001,0.001):
                #Calculate position & velocity at end of hexapole
                thex = d / velocity['vx'] # thex = time spent in traveling to collision
                xd.append(d+d0)
                yd.append((velyf * thex) + ry0)
                zd.append((velzf * thex) + rz0)
                if p.get('pin_pos') != None:
                    if p['pin_pos'] == round(d + 0.001,3):
                        if check_radial_position(p['pin_r'],yd[-1],zd[-1]) == 1: 
                            crashed = True
                            break

            velocity['vy'] = velyf
            velocity['vz'] = velzf

        if check_radial_position(0.00005,yd[-1],zd[-1]) == 1: crashed = True
        
        d2.draw_molecule(xd, yd, zd)
        d3.draw_molecule(xd, yd, zd)
        if crashed != True: count+=1

    print(f'{round((count/(est-skimmer))*100,2)} % of molecules that pass the skimmer made it to collision area.')
    
    #Check if user wants to show plots
    if md['plot_2d_t'] != True: d2.destroy()
    if md['plot_3d_t'] != True: d3.destroy()
    if md['plot_2d_t'] == True or md['plot_3d_t'] == True: drawing.draw_plot()

