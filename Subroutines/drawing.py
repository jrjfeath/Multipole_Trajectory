import matplotlib.pyplot as plt
import numpy as np

class d2_drawing():
    '''
    Make a class to handle drawing the 2d trajectory plots.
    '''
    def __init__(self,scale=1000,colour='black'):
        self.fig = plt.figure()
        self.fig.set_size_inches(10.5, 5.5)
        self.y_plot = self.fig.add_subplot(121)
        self.y_plot.set_xlabel('Distance travelled (mm)')
        self.y_plot.set_ylabel('y-axis travel (mm)')
        self.z_plot = self.fig.add_subplot(122)
        self.z_plot.set_xlabel('Distance travelled (mm)')
        self.z_plot.set_ylabel('z-axis travel (mm)')
        #What are we converting to? m -> mm * 1000, m -> cm *100, etc
        self.scale = scale
        self.colour = colour

    def draw_multipole(self,x,y):
        '''
        Draw in a multipole given x and y coordinates.
        '''
        x *= self.scale
        y *= self.scale 

        self.y_plot.plot(x,y,color=self.colour)
        self.y_plot.plot(x,-y,color=self.colour)
        self.z_plot.plot(x,y,color=self.colour)
        self.z_plot.plot(x,-y,color=self.colour)

    def draw_molecule(self,xd,yd,zd):
        '''
        Draw in the molecule position.
        '''
        #Convert the measurements from meters to mm
        xd = np.array(xd) * self.scale
        yd = np.array(yd) * self.scale
        zd = np.array(zd) * self.scale
        self.y_plot.plot(xd,yd,color=self.colour, alpha=0.2)
        self.z_plot.plot(xd,zd,color=self.colour, alpha=0.2)

    def draw_collision(self,length,d=0.001):
        '''
        Draw the area representing the collision area.
        '''
        #Convert the measurements from meters to mm
        x = np.array([length,length]) * self.scale
        y = np.array([-d, d]) * self.scale
        self.y_plot.plot(x,y,color=self.colour)
        self.z_plot.plot(x,y,color=self.colour)

    def draw_pinhole(self,length,radius):
        '''
        Draw the lines representing the skimmer or pinhole walls.
        '''
        x = np.array([length,length]) * self.scale
        y = np.array([radius, 0.02]) * self.scale
        self.y_plot.plot(x,y,color=self.colour)
        self.y_plot.plot(x,-y,color=self.colour)
        self.z_plot.plot(x,y,color=self.colour)
        self.z_plot.plot(x,-y,color=self.colour)

    def find_poles_and_collision(self,md):
        '''
        Loop through the parameter data to find and draw multipoles and 
        collision centre.
        '''
        cxl = md['lskimmer']
        self.draw_pinhole(cxl,md['skmr_radius'])
        cxl += md['lsource']
        mr = 0
        for m in md['multipole']:
            p = md['multipole'][m]
            #Arrays for our multipole dimensions
            mx = np.arange(0,p['lpole'],0.001) + cxl
            my = np.full((len(mx), 1), p['r0'])
            self.draw_multipole(mx,my)
            if mr < p['r0']: mr = p['r0']
            #Add multipole distances to total
            cxl+=p['lpole']
            if p.get('pin_pos') != None:
                self.draw_pinhole(cxl + p['pin_pos'],p['pin_r'])
            cxl+=p['ldist']
        #Set bounds for drawing
        mr = (mr * 1000) + 0.5
        self.y_plot.set_ylim(ymin=-mr, ymax=mr)
        self.z_plot.set_ylim(ymin=-mr, ymax=mr)
        #Add collision distance to total
        cxl+=md['lcollision']
        self.draw_collision(cxl)

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

    def multipole_radius(self,md):
        '''
        Determine the largest multipole radius and limit the axes.
        '''
        rm = 0
        for m in md['multipole']:
            p = md['multipole'][m]
            if rm < p['r0']: rm = p['r0']
        rm *= self.scale
        self.ax3d.set_ylim(-rm,rm)
        self.ax3d.set_zlim(-rm,rm)

    def draw_molecule(self,xd,yd,zd):
        #Convert the measurements from meters to mm
        xd = np.array(xd) * self.scale
        yd = np.array(yd) * self.scale
        zd = np.array(zd) * self.scale
        self.ax3d.plot3D(xd, yd, zd, c=self.colour, alpha=0.2)

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