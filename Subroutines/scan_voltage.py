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

    #Loop through a specified number of molecules
    for _ in range(d['vest']):
        #Density index, used for tracking which index to write to
        di = 0
        #Calculate an initial trajectory for each of our molecule
        velocity = skimmer_trajectory(d)
        #Make an empty array to hold velocity data after it is modified
        vra = np.zeros((len(d['V']),4))
        #Make an array to track successful passes
        vsa = np.zeros((len(d['V']),1))
        #fetch the initial distance from skimmer to multipole
        lsource = d['lsource']
        #Loop through each multipole
        for mindex, m in enumerate(d['multipole']):
            #fetch multipole parameters
            p = d['multipole'][m]
            sqc , invc, r0, ldist = p['sqc'], p['invc'],  p['r0'], p['ldist']
            #Loop through each voltage
            for vindex, voltage in enumerate(d['V']):
                #Track how much di should increase by this loop
                dii = 0

                #On the first multipole we want to use our initial velocities
                if mindex == 0:
                    #Calculate our trajectory at the end of the multipole
                    info = multipole_trajectory(
                        lsource,
                        p['lpole'],
                        velocity['vx'],
                        velocity['vy'],
                        velocity['vz'],
                        velocity['ry'],
                        velocity['rz'],
                        sqc[vindex],
                        invc[vindex],
                        r0
                    )
                
                #After the first multipole we want to use the values the values
                #calculated from the multipole prior to this one
                else:
                    #If the previous loop failed to exit the multipole, skip
                    if vsa[vindex] == 0: continue
                    #Calculate our trajectory at the end of the multipole
                    info = multipole_trajectory(
                        lsource,
                        p['lpole'],
                        velocity['vx'],
                        vra[vindex][0],
                        vra[vindex][1],
                        vra[vindex][2],
                        vra[vindex][3],
                        sqc[vindex],
                        invc[vindex],
                        r0
                    )
                
                #Write entry data to densities
                densities[vindex][di+dii] += info['enter']
                dii+=1

                #overwrite velocity and position data
                vra[vindex] = [info['velyf'],info['velzf'],info['ryhex'],info['rzhex']]
                vsa[vindex] = info['exit']

                #If molecule/ion does not exit multipole continue to next voltage
                if info['exit'] == 0: continue

                #If the multipole has a pinhole after it calculate success
                #Overwrite the exit parameter to make data storage less cumbersome
                if p.get('pin_pos') != None:
                    #If the molecule/ion passes region returns 1
                    info['exit'] = check_collision(
                        velocity['vx'],
                        vra[vindex][0], 
                        vra[vindex][1], 
                        vra[vindex][2], 
                        vra[vindex][3],
                        p['pin_pos'],
                        p['pin_r']
                    )

                #overwrite the position data on if it travels through pinhole
                vsa[vindex] = info['exit']

                #Write exit data to densities
                densities[vindex][di+dii] += info['exit']
                dii+=1

                #If we are on the last multipole calculate final position at
                #the collision centre if molecule exited multipole
                if mindex == len(d['multipole']) - 1 and info['exit'] == 1:
                    densities[vindex][-1] += check_collision(
                        velocity['vx'],
                        vra[vindex][0], 
                        vra[vindex][1], 
                        vra[vindex][2], 
                        vra[vindex][3],
                        d['lcollision']
                    )

            #Increase the density index, either by 2 or 3
            di += dii

            #After we finish calculating trajectories change the distance
            #to the next part of the setup
            lsource = ldist

    for index, row in enumerate(densities):
        #print(d['V'][index],row)
        pass

    #Plot the data
    x = [x for x in d['V']]
    y = (densities[::,-1] / d['vest']) * 100
    plot.draw(x,y)

    #Get the maximum transmission and append to plot
    maximum = max(densities[::,-1])
    mv = (maximum / d['vest']) * 100
    index = np.where(densities[::,-1]==maximum)[0][-1]
    label = f'{x[index]}, {round(mv,2)}%'
    plot.annotate(label,x[index],mv)
    
    drawing.draw_plot()
    #Destroy the plot in the event user runs multiple scan commands
    plot.destroy()
