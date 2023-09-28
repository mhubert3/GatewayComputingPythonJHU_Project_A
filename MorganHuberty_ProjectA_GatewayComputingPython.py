#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 19:59:28 2023

@author: morganhuberty
"""

# simVis1 file imported below

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 13:43:42 2023

@author: siamak
"""



import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
#import time



class visualize:
    
    np.random.seed(123456789)
    
    colorDict = {0: 'green', 1:'red', 2:'gray', 3:'black'}
    faceColorDict = {0: (0, 1, 0, .3) , 1: (1, 0, 0, .3), \
                     2: (0.5, 0.5, 0.5, .3),3: (0, 0, 0, .8)}
    
    
    def __init__(self,figwWidth=5,figHeight=5,pauseTime=0.001, 
                 Ax1xmin=0,Ax1xmax = 1, Ax1ymin = 0, Ax1ymax = 1):
        
        
        self.figwWidth = figwWidth
        self.figHeight = figHeight
        
        self.Ax1xmin = Ax1xmin
        self.Ax1xmax = Ax1xmax
        self.Ax1ymin = Ax1ymin
        self.Ax1ymax = Ax1ymax
        
        self.pauseTime = pauseTime
        
        self.f, self.ax1 = plt.subplots(1, 1,figsize=(self.figwWidth,self.figHeight ))
        self.setAx1Limits()
        
        
        
    def circle(self,x, y, radius,colorCode):
        """
        Draw one circle with center placed at x and y and radius value of radius. Sets 
        edge and facecollors based on the colorCode.

        Parameters
        ----------
        x : float
            x-coordinate of the center
        y : float
            y-coordinate of the center
        radius : float
            radius of the circle
        colorCode : integer
           defines the face and edge colors  using class dictionary  attributes.

        Returns
        -------
        None.

        """
        
    
        circle = Circle((x, y), radius , edgecolor=self.colorDict[colorCode] ,facecolor=self.faceColorDict[colorCode])
        self.ax1.add_artist(circle)
        self.setAx1Limits()
        

            
            
    def getXminMax(self):
        """
        gets axis 1 x-min and x-max values

        Returns
        -------
        None.

        """
        
        self.Ax1xmin, self.Ax1xmax = self.ax1.axes.get_xlim() 


    def getYminMax(self):
        """
        gets axis 1 y-min and y-max values

        Returns
        -------
        None.

        """
        
        self.Ax1ymin, self.Ax1ymax = self.ax1.axes.get_xlim() 
        
    def setAx1Limits(self):
        """
        sets axis 1 x and y limits

        Returns
        -------
        None.

        """
        
        self.ax1.set_xlim([self.Ax1xmin, self.Ax1xmax])
        self.ax1.set_ylim([self.Ax1ymin, self.Ax1ymax])
        
    
 
    
    
    def plotPause(self):
        """
        Pauses the plot for the duration of pauseTime

        Returns
        -------
        None.

        """
        
        plt.pause(self.pauseTime)
        
    def axis1Clear(self):
        """
        Clears axis 1 in subplot (left panel)

        Returns
        -------
        None.

        """
        
        self.ax1.cla()
        
    # For debugging
    def addText(self,circlesX,circlesY):
        
        font = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 8
        }

        # #add text with custom font
        # for i, j, k in zip(Agents.x,Agents.y,np.arange(0,Agents.getNumberOfCircles()+1)):
        #     self.ax1.text(i, j, k, fontdict=font)
            
        #add text with custom font
        for i, j, k in zip(circlesX,circlesY,np.arange(0,2+1)):
            self.ax1.text(i , j - 0.03, k+1, fontdict=font)
            
    def addLine(self,circlesX,circlesY,line_style = '--',marker='.'):
        
        self.ax1.plot(circlesX,circlesY,linestyle = line_style, marker = marker)
        self.setAx1Limits()
        
    def drawVector(self,initial_point,vector):
        
        self.ax1.arrow(initial_point[0],initial_point[1], 
                        vector[0],vector[1],
                        head_width=0.01,head_length=0.01, color='r')
        #angles='xy', scale_units='xy', 
        self.setAx1Limits()
            
            
if __name__ == '__main__':
    
    pass
    
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 16:05:56 2023


"""

# projectA file imported below

import math




def main():
    
    # initialize parameters
    
    radius =0.15
    dt =.2
    
    # simulation steps
    simSteps = 800
    

    # Default is Ax1xmin= 0,Ax1xmax = 1, Ax1ymin = 0, Ax1ymax = 1
    vis = visualize() 
    

    x1 = 0.2
    y1 = 0.25
    
    v1x = -0.05
    v1y = -0.05
    
    x2 = 0.7
    y2 = 0.7
    
    
    v2x = 0.05
    v2y = 0.05
    
    
    # run simulation
    simulate(simSteps,vis,x1,y1,x2,y2,v1x,v1y,v2x,v2y,dt,radius)
    


# DO NOT CHANGE THIS FUNCTION
def boundary_locations(vis,radius):
    
    """
    Calculate the boundary locations within the visualization window.
    This function calculates the boundary coordinates for a containing box 
    within the visualization window to ensure that particles 
    do not exceed these boundaries.

    Parameters
    ----------
    vis : visualization object
        The visualization object used for determining window boundaries.
    radius : float
        The radius of the particles.

    Returns
    -------
    xLow : float
        The lower bound of the containing box's threshold to hit the vertical left wall.
    yLow : float
        The lower bound of the containing box's threshold to hit the horizontal bottom wall.
    xHigh : float
        The upper bound of the containing box's threshold to hit the vertical right wall.
    yHigh : float
        The upper bound of the containing box's threshold to hit the horizontal top wall.

    """
    
    # Adjust the boundary according to the 
    # circle radius and display window size
    # Assume circles have the same radius.
    xLow = vis.Ax1xmin + radius 
    xHigh = vis.Ax1xmax - radius 
    
    yLow = vis.Ax1ymin + radius 
    yHigh = vis.Ax1ymax - radius 
    
    return xLow, yLow, xHigh, yHigh
    
    


def updateX(x,vx,dt):
    
    """
    Update the x-coordinate of the center of a circle's position
    based on its velocity over a specified time step.

    Parameters
    ----------
    x : float
        The current x-coordinate of a circle.
    vx : float
        The current x-component velocity of a circle.
    dt : float
        Time step size for the simulation.

    Returns
    -------
    xnew : float
        Updated x-coordinate of a circle.

    """
    
    # Update the x-coordinate based on the velocity and time step
    xnew = x + (vx * dt)

    # Round the result to 3 decimal points
    xnew = round(xnew, 3)

    return xnew

def updateY(y,vy,dt):
    
    """
    Update the y-coordinate of the center of a circle's position
    based on its velocity over a specified time step.

    Parameters
    ----------
    y : float
        The current y-coordinate of a circle.
    vy : float
        The current y-component velocity of a circle.
    dt : float
        Time step size for the simulation.

    Returns
    -------
    ynew : float
        Updated y-coordinate of a circle.

    """
    
    # Update the y-coordinate based on the velocity and time step
    ynew = y + (vy * dt)

    # Round the result to 3 decimal points
    ynew = round(ynew, 3)

    return ynew

def boxCollision(x,y,vx,vy,xLow,xHigh,yLow,yHigh):
    
    """
    Determine whether circle 1 or circle 2 is in contact with a threshold of 
    the containing box.

    Parameters
    ----------
    x : float
        The current x-coordinate of a circle.
    y : float
        The current y-coordinate of a circle.
    vx : float
        The current x-component velocity of a circle.
    vy : float
        The current y-component velocity of a circle.
    xLow : float
        An x-component threshold defining the containing box.
    xHigh : float
        An x-component threshold defining the containing box.
    yLow : float
        An y-component threshold defining the containing box.
    yHigh : float
        An y-component threshold defining the containing box.

    Returns
    -------
    Other returns are defined in parameters section.

    """
    
    # Check if the circle hits the left wall
    if x - vx < xLow:
        x = xLow  # Set x to the left wall
        vx = -vx  # Reflect velocity in the x-direction

    # Check if the circle hits the right wall
    if x + vx > xHigh:
        x = xHigh  # Set x to the right wall
        vx = -vx  # Reflect velocity in the x-direction

    # Check if the circle hits the bottom wall
    if y - vy < yLow:
        y = yLow  # Set y to the bottom wall
        vy = -vy  # Reflect velocity in the y-direction

    # Check if the circle hits the top wall
    if y + vy > yHigh:
        y = yHigh  # Set y to the top wall
        vy = -vy  # Reflect velocity in the y-direction

    return x, y, vx, vy
        
def Overlap(x1,y1,radius1, x2,y2,radius2):
    
    """
    Determine whether circle 1 and circle 2 overlap or collide.

    Parameters
    ----------
    x1 : float
        Initial x-coordinate of the center of circle 1.
    y1 : float
        Initial y-coordinate of the center of circle 1.
    radius1 : float
        Radius of circle 1.
    v1x : float
        Initial x-component of the velocity of circle 1.
    v1y : float
        Initial y-component of the velocity of circle 1.
    x2 : float
        Initial x-coordinate of the center of circle 2.
    y2 : float
        Initial y-coordinate of the center of circle 2.
    radius2 : float
        Radius of circle 2.
    v2x : float
        Initial x-component of the velocity of circle 2.
    v2y : float
        Initial y-component of the velocity of circle 2.

    Returns
    -------
    x1new : float
        Post-collision x-coordinate of the center of circle 1.
    y1new : float
        Post-collision y-coordinate of the center of circle 1.
    v1xnew : float
        Post-collision x-component of the velocity of circle 1.
    v1ynew : float
        Post-collision y-component of the velocity of circle 1.
    x2new : float
        Post-collision x-coordinate of the center of circle 2.
    y2new : float
        Post-collision y-coordinate of the center of circle 2.
    v2xnew : float
        Post-collision x-component of the velocity of circle 2.
    v2ynew : float
        Post-collision y-component of the velocity of circle 2.

    """
    
    # Calculate the distance between the centers of the two circles
    distance = get_distance(x1, y1, x2, y2)

    # Check if the circles overlap
    if distance < (radius1 + radius2):
        return True
    else:
        return False

def get_unit_direction(x1,y1,x2,y2):
    
    """
    Calculates the unit direction vector between two points represented by
    coordinates of the centers or circles 1 and 2.

    Parameters
    ----------
    x1 : float
        Initial x-coordinate of the center of circle 1.
    y1 : float
        Initial y-coordinate of the center of circle 1.
    x2 : float
        Initial x-coordinate of the center of circle 2.
    y2 : float
        Initial y-coordinate of the center of circle 2.

    Returns
    -------
    x-unit : float
        X-component of unit direction vecotor from centers of circle 1 to circle 2.
    Y-unit : float
        Y-component of unit direction vecotor from centers of circle 1 to circle 2.

    """
    
    # Calculate the denominator using the get_distance function
    denominator = math.sqrt(((x1 - x2) ** 2) + ((y1 - y2) ** 2))

    # Check if the two points coincide
    if math.isclose(denominator, 0, rel_tol=1e-9):
        return math.nan, math.nan

    # Calculate the unit direction vector components
    x_unit = (x1 - x2) / denominator
    y_unit = (y1 - y2) / denominator

    # Round the values to 3 decimal points
    x_unit = round(x_unit, 3)
    y_unit = round(y_unit, 3)

    print(f"The unit direction vector from B to A is ({x_unit}, {y_unit})")

    return x_unit, y_unit

def get_distance(x1,y1,x2,y2):
    
    """
    Calculates the Euclidean distance between two points represented by
    coordinates of the centers or circles 1 and 2.

    Parameters
    ----------
    x1 : float
        Initial x-coordinate of the center of circle 1.
    y1 : float
        Initial y-coordinate of the center of circle 1.
    x2 : float
        Initial x-coordinate of the center of circle 2.
    y2 : float
        Initial y-coordinate of the center of circle 2.

    Returns
    -------
    distance : float
        Calculation of distance between centers of circle 1 and circle 2.

    """
    
    # Calculate the squared differences
    dx = (x2 - x1) ** 2
    dy = (y2 - y1) ** 2
    
    # Calculate the Euclidean distance
    distance = math.sqrt(dx + dy)
    
    print(f"The Euclidean distance between A({x1}, {y1}) and B({x2}, {y2}) is {distance}")
    
    # Round the distance to 3 decimal points
    return round(distance, 3)
    
def dot_product(x1,y1,x2,y2):
    
    """
    Calculates the dot product of two 2D vectors represented by
    coordinates of the centers or circles 1 and 2.

    Parameters
    ----------
    x1 : float
        Initial x-coordinate of the center of circle 1.
    y1 : float
        Initial y-coordinate of the center of circle 1.
    x2 : float
        Initial x-coordinate of the center of circle 2.
    y2 : float
        Initial y-coordinate of the center of circle 2.

    Returns
    -------
    result : float
        Calculation of dot product.

    """
    
    # Calculate the dot product
    result = (x1 * x2) + (y1 * y2)
    
    print(f"The dot product of v1({x1}, {y1}) and v2({x2}, {y2}) is {result}")
    
    return result

def update_collision_velocity(x1,y1,u1x,u1y, x2,y2,u2x,u2y):
    
    """
    Calculate the updated velocities of the circles after each collison.

    Parameters
    ----------
    x1 : float
        Initial x-coordinate of the center of circle 1.
    y1 : float
        Initial y-coordinate of the center of circle 1.
    u1x : float
        X-component of the current velocity vector of circle 1.
    u1y : float
        Y-component of the current velocity vector of circle 1.
    x2 : float
        Initial x-coordinate of the center of circle 2.
    y2 : float
        Initial y-coordinate of the center of circle 2.
    u2x : float
        X-component of the current velocity vector of circle 2.
    u2y : float
        Y-component of the current velocity vector of circle 2.

    Returns
    -------
    V1x : float
        Post-collision updated x-component velocity of circle 1.
    V1y : float
        Post-collision updated y-component velocity of circle 1.
    V2x : float
        Post-collision updated x-component velocity of circle 2.
    V2y : float
        Post-collision updated y-component velocity of circle 2.

    """
    
    # Calculate the Euclidean distance between the centers of the circles
    distance = math.sqrt((x2 - x1) ** 2) + ((y2 - y1) ** 2)
    
    # Check if the circles overlap (collision)
    if distance == 0:
        return -u1x, -u1y, -u2x, -u2y  # Velocities are reversed

    # Calculate the updated velocities based on the collision equations
    V1x = u1x + (-1 * ((u1x - u2x) * (x1 - x2) + (u1y - u2y) * (y1 - y2)) / (distance ** 2)) * (x1 - x2)
    V1y = u1y + (-1 * ((u1x - u2x) * (x1 - x2) + (u1y - u2y) * (y1 - y2)) / (distance ** 2)) * (y1 - y2)
    V2x = u2x + (-1 * ((u2x - u1x) * (x2 - x1) + (u2y - u1y) * (y2 - y1)) / (distance ** 2)) * (x2 - x1)
    V2y = u2y + (-1 * ((u2x - u1x) * (x2 - x1) + (u2y - u1y) * (y2 - y1)) / (distance ** 2)) * (y2 - y1)

    # Round the updated velocities to 3 decimal points
    V1x = round(V1x, 3)
    V1y = round(V1y, 3)
    V2x = round(V2x, 3)
    V2y = round(V2y, 3)

    return V1x, V1y, V2x, V2y

def circleCollision(x1,y1,radius1,v1x,v1y, x2,y2,radius2,v2x,v2y):
   
    """
    Calculate the positions and velocities of the circles after each collison.

    Parameters
    ----------
    x1 : float
        Initial x-coordinate of the center of circle 1.
    y1 : float
        Initial y-coordinate of the center of circle 1.
    radius1 : float
        Radius of circle 1.
    v1x : float
        Initial x-component of the velocity of circle 1.
    v1y : float
        Initial y-component of the velocity of circle 1.
    x2 : float
        Initial x-coordinate of the center of circle 2.
    y2 : float
        Initial y-coordinate of the center of circle 2.
    radius2 : float
        Radius of circle 2.
    v2x : float
        Initial x-component of the velocity of circle 2.
    v2y : float
        Initial y-component of the velocity of circle 2.

    Returns
    -------
    x1new : float
        Post-collision x-coordinate of the center of circle 1.
    y1new : float
        Post-collision y-coordinate of the center of circle 1.
    v1xnew : float
        Post-collision x-component of the velocity of circle 1.
    v1ynew : float
        Post-collision y-component of the velocity of circle 1.
    x2new : float
        Post-collision x-coordinate of the center of circle 2.
    y2new : float
        Post-collision y-coordinate of the center of circle 2.
    v2xnew : float
        Post-collision x-component of the velocity of circle 2.
    v2ynew : float
        Post-collision y-component of the velocity of circle 2.
    Other returns are defined in parameters section.

    """
    
    # Calculate the distance between the centers of the circles
    distance = math.sqrt(((x2 - x1) ** 2) + ((y2 - y1) ** 2))
    
    # Check if the circles overlap (collision)
    if distance < radius1 + radius2:
        # Calculate the unit direction vector along the centers of the circles
        x_unit, y_unit = get_unit_direction(x1, y1, x2, y2)
        
        # Calculate the displacement needed to separate the circles
        displacement = (radius1 + radius2 - distance)
        
        # Update the positions of the circles
        x1New = x1 + displacement * x_unit
        y1New = y1 + displacement * y_unit
        x2New = x2 - displacement * x_unit
        y2New = y2 - displacement * y_unit
        
        # Calculate the updated velocities using the reflection formula
        dot_product = (v1x - v2x) * (x1 - x2) + (v1y - v2y) * (y1 - y2)
        
        v1xNew = v1x - (2 * dot_product / (radius1 ** 2)) * (x1 - x2)
        v1yNew = v1y - (2 * dot_product / (radius1 ** 2)) * (y1 - y2)
        v2xNew = v2x - (2 * dot_product / (radius2 ** 2)) * (x2 - x1)
        v2yNew = v2y - (2 * dot_product / (radius2 ** 2)) * (y2 - y1)
        
        # Round the values to 3 decimal points
        x1New = round(x1New, 3)
        y1New = round(y1New, 3)
        x2New = round(x2New, 3)
        y2New = round(y2New, 3)
        v1xNew = round(v1xNew, 3)
        v1yNew = round(v1yNew, 3)
        v2xNew = round(v2xNew, 3)
        v2yNew = round(v2yNew, 3)
        
        print("New Positions and Velocities After Collision:")
        print("Circle 1:")
        print(f"Position: ({x1New:.2f}, {y1New:.2f})")
        print(f"Velocity: ({v1xNew:.2f}, {v1yNew:.2f})")
        print("Circle 2:")
        print(f"Position: ({x2New:.2f}, {y2New:.2f})")
        print(f"Velocity: ({v2xNew:.2f}, {v2yNew:.2f})")
        
        return x1New, y1New, v1xNew, v1yNew, x2New, y2New, v2xNew, v2yNew
   
    else:
        # No collision, return the original values
        return x1, y1, v1x, v1y, x2, y2, v2x, v2y



# DON NOT CHANGE THIS FUNCTION
def simulate(simSteps,vis,x1,y1,x2,y2,v1x,v1y,v2x,v2y,dt,radius):
    
    """
    Run a particle simulation for a specified number of time steps.

    Parameters
    ----------
    simSteps : int
        Number of time steps to run the simulation.
    vis : visualize
        An instance of the 'visualize' class for visualization.
    x1 : float
        Initial x-coordinate of the center of circle 1.
    y1 : float
        Initial y-coordinate of the center of circle 1.
    x2 : float
        Initial x-coordinate of the center of circle 2.
    y2 : float
        Initial y-coordinate of the center of circle 2.
    v1x : float
        Initial x-component of the velocity of circle 1.
    v1y : float
        Initial y-component of the velocity of circle 1.
    v2x : float
        Initial x-component of the velocity of circle 2.
    v2y : float
        Initial y-component of the velocity of circle 2.
    dt : float
        Time step size for the simulation.
    radius : float
        Radius of both particles.

    Returns
    -------
    None

    """
    
    # This function runs a two-particle simulation for a specified number of 
    # time steps. It updates the positions of two circles and checks for 
    # their collisions with both the container boundary and each other. 
    # The simulation is visualized using the'visualize' class instance 
    # provided as 'vis'.
    
    
    # Identify the boundary location for the containing window
    # Assumption: both particles have the same radius.
    
    xLow,yLow,xHigh, yHigh = boundary_locations(vis,radius) 
    
    for i in range(simSteps):
        
        # update x, y position of particles
        x1 = updateX(x1,v1x,dt)
        y1 = updateY(y1,v1y,dt)
        
        x2 = updateX(x2,v2x,dt)
        y2 = updateY(y2,v2y,dt)
            
        
                    
        # Check for box and circle collisions
        
        x1,y1,v1x,v1y = boxCollision(x1,y1,v1x,v1y,xLow,xHigh,yLow,yHigh)
        
        
        
        x2,y2,v2x,v2y = boxCollision(x2,y2,v2x,v2y,xLow,xHigh,yLow,yHigh)
        
        
        x1,y1,v1x,v1y, x2,y2,v2x,v2y = circleCollision\
                                      (x1,y1,radius,v1x, \
                                        v1y, x2,y2,radius,v2x,v2y)
                    
        
                    
        # draw circles
        vis.circle(x1, y1, radius, 0)
        vis.circle(x2, y2, radius, 1)
        
                    
        # pause plots and clear window axis 1
        vis.plotPause()
        vis.axis1Clear()
    
    # redraw circles after last iteration
    vis.circle(x1, y1, radius, 0)
    vis.circle(x2, y2, radius, 1)
        
if __name__ == '__main__': 
    
    # Call the main function to excute the simulation
    # For testing purpose comment main() and call another 
    # function.
    main()
    
    
    

    