import matplotlib.pyplot as plt
import numpy as np

class hist_coll():
    def __init__(self,scale=1000,colour='black',mc='purple'):
        self.fig = plt.figure()
        self.fig.set_size_inches(5.5, 5.5)
        self.rad = self.fig.add_subplot(111)
        self.scale = scale

    def draw_hist(self,x,y):
        self.rad.plot(x * self.scale,y)
        self.rad.set_xlabel('Radius (mm)')
        self.rad.set_ylabel('Intensity')
            
class d2_drawing():
    '''
    Make a class to handle drawing the 2d trajectory plots.
    '''
    def __init__(self,scale=1000,colour='black',mc='purple'):
        self.fig = plt.figure()
        self.fig.set_size_inches(5.5, 5.5)
        self.y_plot = self.fig.add_subplot(111)
        self.y_plot.set_xlabel('Distance travelled (mm)')
        self.y_plot.set_ylabel('Radius (mm)')
        #What are we converting to? m -> mm * 1000, m -> cm *100, etc
        self.scale = scale
        self.colour = colour
        self.molecule_color = mc
        self.length = 0

    def draw_multipole(self,x,y,label):
        '''
        Draw in a multipole given x and y coordinates.
        '''
        x *= self.scale
        y *= self.scale 

        self.y_plot.plot(x,y,color=self.colour)
        self.y_plot.text(x[0],y[0]+0.1,label)
        self.y_plot.plot(x,-y,color=self.colour)

    def draw_molecule(self,xd,yd):
        '''
        Draw in the molecule position.
        '''
        #Convert the measurements from meters to mm
        xd = np.array(xd) * self.scale
        yd = np.array(yd) * self.scale
        self.y_plot.plot(xd,yd,color=self.molecule_color, alpha=0.02)

    def draw_collision(self,length,d=0.001):
        '''
        Draw the area representing the collision area.
        '''
        #Convert the measurements from meters to mm
        x = np.array([length,length]) * self.scale
        y = np.array([-d, d]) * self.scale
        self.y_plot.plot(x,y,color=self.colour)
        self.y_plot.text(x[0] + 15,y[0] - 0.5,'Collision Region',rotation=90)

    def draw_pinhole(self,length,radius,label):
        '''
        Draw the lines representing the skimmer or pinhole walls.
        '''
        x = np.array([length,length]) * self.scale
        y = np.array([radius, 0.02]) * self.scale
        self.y_plot.plot(x,y,color=self.colour)
        self.y_plot.text(x[0]+1,y[0]+1,label)
        self.y_plot.plot(x,-y,color=self.colour)

    def find_poles_and_collision(self,md):
        '''
        Loop through the parameter data to find and draw multipoles and 
        collision centre.
        '''
        cxl = md['lskimmer']
        self.draw_pinhole(cxl,md['skmr_radius'],'Skimmer')
        cxl += md['lsource']
        mr = 0
        for m in md['multipole']:
            p = md['multipole'][m]
            #Arrays for our multipole dimensions
            mx = np.arange(0,p['lpole'],0.001) + cxl
            my = np.full((len(mx), 1), p['r0'])
            self.draw_multipole(mx,my,m)
            if mr < p['r0']: mr = p['r0']
            #Add multipole distances to total
            cxl+=p['lpole']
            if p.get('pin_pos') != None:
                self.draw_pinhole(cxl + p['pin_pos'],p['pin_r'],'Pinhole')
            cxl+=p['ldist']
        #Set bounds for drawing
        mr = (mr * 1000) + 0.5
        self.y_plot.set_ylim(ymin=-mr, ymax=mr)
        #Add collision distance to total
        cxl+=md['lcollision']
        self.draw_collision(cxl,md['crr'])
        self.length = cxl

    def destroy(self):
        plt.close(self.fig)

class d3_drawing():
    '''
    Make a class to handle drawing the 3d trajectory plots.
    '''
    def __init__(self,scale=1000,colour='black'):
        self.fig = plt.figure()
        self.ax3d = plt.axes(projection='3d')
        self.ax3d.set_xlabel('Distance travelled (mm)')
        self.ax3d.set_ylabel('y-axis travel (mm)')
        self.ax3d.set_zlabel('z-axis travel (mm)')
        self.scale = scale
        self.colour = colour

    def draw_multipole(self,center_x, center_y, radius, height_z):
        z = np.linspace(0, height_z, 50)
        theta = np.linspace(0, 2*np.pi, 50)
        theta_grid, z_grid = np.meshgrid(theta, z)
        x_grid = radius * np.cos(theta_grid) + center_x
        y_grid = radius * np.sin(theta_grid) + center_y
        return x_grid,y_grid,z_grid

    def multipole_radius(self,md):
        '''
        Determine the largest multipole radius and limit the axes.
        '''
        rm = 0
        for m in md['multipole']:
            p = md['multipole'][m]
            if rm < p['r0']: rm = p['r0']
            r1 = round(p['d0'] * 0.5628,1) #Radius of multipole rod
            r0 = p['r0'] * 1000
            length = p['lpole'] * 1000
            degree = 360 / p['size']
            centre = round(r0 + r1/2,1)
            for i in range(p['size']):
                x_centre = round(centre * np.cos(np.radians(degree * i)),1)
                y_centre = round(centre * np.sin(np.radians(degree * i)),1)
                Xc,Yc,Zc = self.draw_multipole(x_centre,y_centre,r1 / 2,length)
                self.ax3d.plot_surface(Zc, Xc, Yc, alpha=0.5, color = 'blue')

    def draw_molecule(self,xd,yd,zd):
        #Convert the measurements from meters to mm
        xd = np.array(xd) * self.scale
        yd = np.array(yd) * self.scale
        zd = np.array(zd) * self.scale
        self.ax3d.plot3D(xd, yd, zd, c=self.colour, alpha=0.05)

    def destroy(self):
        plt.close(self.fig)

class voltage_plot():
    '''
    Make a class to handle drawing the voltage plot.
    '''
    def __init__(self,colour='black'):
        title = 'What percent of molecules are transmitted to collision region?'

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlabel('Voltage (V)')
        self.ax.set_ylabel('Percent Transmission (%)')
        self.ax.set_title(title)
        self.ax.set_ylim([0,101])
        self.colour = colour

    def draw(self,x,y):
        self.ax.plot(x,y,self.colour)

    def annotate(self,label,x,y):
        self.ax.annotate(label, (x,y))

    def destroy(self):
        plt.close(self.fig)

def draw_plot():
    '''Generic function for showing plots.'''
    plt.show()