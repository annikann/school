import sys
sys.path.append('/Users/annikacarlson/Documents/ControlTheory/smdSim')
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
import numpy as np 
import parameters.smdParam as P

class smdAnimation:
    def __init__(self):
        self.flag_init = True  # Used to indicate initialization
        # Initialize a figure and axes object
        self.fig, self.ax = plt.subplots()
        # Initializes a list of objects (patches and lines)
        self.handle = []
        # Specify the x,y axis limits
        plt.axis([-2, 2, -0.1, 3])
        # Draw line for the ground
        plt.plot([-3, 3], [0, 0], 'b--')
        # label axes
        plt.xlabel('z')

    def update(self, state):
        z = state[0][0]  # Horizontal position of cart, m
        # draw plot elements
        self.draw_cart(z)
        self.ax.axis('equal')
        # Set initialization flag to False after first call
        if self.flag_init == True:
            self.flag_init = False

    def draw_cart(self, z):
        # specify bottom left corner of rectangle
        x = z-P.w/2.0
        y = P.gap
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