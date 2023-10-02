import sys
sys.path.append('/Users/annikacarlson/Documents/school/controltheory/library')
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
import numpy as np 
import library.vtolParam as P

class vtolAnimation:
    def __init__(self):
        self.flag_init = True  # Used to indicate initialization
        # Initialize a figure and axes object
        self.fig, self.ax = plt.subplots()
        # Initializes a list of objects (patches and lines)
        self.handle = []
        # Specify the x,y axis limits
        plt.axis([-10, 10, -1, 20])
        # Draw line for the ground
        plt.plot([-10, 10], [0, 0], 'b--')
        # label axes
        plt.xlabel('z')

    def update(self, state):
        z = state[0][0]      # Horizontal position of vtol, m
        h = state[2][0]      # Vertical position of vtol, m
        theta = state[4][0]  # Angle of vtol, rads
        # draw plot elements
        self.draw_body(z, h)
        self.draw_right(z, h, theta)
        self.draw_left(z, h, theta)
        # self.ax.axis('equal')
        # Set initialization flag to False after first call
        if self.flag_init == True:
            self.flag_init = False

    def draw_body(self, z, h):
        # specify bottom left corner of rectangle
        x = z - P.w/2.0
        y = h - P.h/2.0
        corner = (x, y)
        # create rectangle on first call, update on subsequent calls
        if self.flag_init is True:
            # Create the Rectangle patch and append its handle
            # to the handle list
            self.handle.append(
                mpatches.Rectangle(corner, P.w, P.h, fc='blue', ec='black'))
            # Add the patch to the axes
            self.ax.add_patch(self.handle[0])
        else:
            self.handle[0].set_xy(corner)  # Update patch

    def draw_right(self, z, h, theta):
        # Define center of prop
        x = [z + ((P.w/2) + (P.d*np.cos(-theta)))]
        y = [h - (P.d*np.sin(-theta))]
        center = (x,y)
        # Create circle on first call, update on subsequent calls
        if self.flag_init is True:
            # Create the circle patch and append its handle to the list
            self.handle.append(
                mpatches.CirclePolygon(center, radius=0.1, resolution=15, fc='blue', ec='black'))
            # Add the path to the axis
            self.ax.add_patch(self.handle[1])
        else: self.handle[1].xy = center

    def draw_left(self, z, h, theta):
        # Define center of prop
        x = [z - ((P.w/2) + (P.d*np.cos(-theta)))]
        y = [h + (P.d*np.sin(-theta))]
        center = (x,y)
        # Create circle on first call, update on subsequent calls
        if self.flag_init is True:
            # Create the circle patch and append its handle to the list
            self.handle.append(
                mpatches.CirclePolygon(center, radius=0.1, resolution=15, fc='blue', ec='black'))
            # Add the path to the axis
            self.ax.add_patch(self.handle[2])
        else: self.handle[2].xy = center
