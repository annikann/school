import sys
sys.path.append('/Users/annikacarlson/Documents/school/controltheory/library')
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
import numpy as np 
import library.bbParam as P

class bbAnimation:
    def __init__(self, limits, multfigs=False):
         # set up plot
        self.fig = plt.figure(1)
        if multfigs == True:
            self.ax = self.fig.add_subplot(1, 2, 1)
        else:
            self.ax = self.fig.add_subplot(1, 1, 1)
        # draw ground
        plt.plot([-limits, limits], [0, 0])
        # for lsit of objects
        self.handle = []
        # init flag
        self.flag_init = True
        # set limits
        self.limits = limits

    def update(self, state):
        z = state[0][0]      # Horizontal position of ball along beam, m
        theta = state[2][0]  # Angle of beam, rads
        
        # draw plot elements
        self.draw_ball(z, theta)
        self.draw_beam(theta)

        # Set axes and limits
        self.ax.set_aspect("equal")
        self.ax.set_ylim(top=self.limits*(5/6), bottom=-self.limits*(5/6))
        self.ax.set_xlim(left=-0.01, right=self.limits)
        # Set initialization flag to False after first call
        if self.flag_init == True:
            self.flag_init = False

    def draw_ball(self, z, theta):
        # Define center of ball
        x = z*np.cos(theta) + P.r*np.cos(theta + np.deg2rad(90))
        y = z*np.sin(theta) + P.r*np.sin(theta + np.deg2rad(90))
        center = (x, y)
        # Create circle on first call, update on subsequent calls
        if self.flag_init is True:
            # Create the circle patch and append its handle to the list
            self.handle.append(mpatches.CirclePolygon(center, radius=P.r, resolution=15, fc='blue', ec='blue'))
            # Add the patch to the axes
            self.ax.add_patch(self.handle[0])
        else:
            self.handle[0].xy = center  # Update patch

    def draw_beam(self, theta):
        # Create rectangle on first call, update on subsequent calls
        if self.flag_init is True:
            # Create the rectangle patch and append its handle to the list
            self.handle.append(mpatches.Rectangle((0,0), P.l, 0.0055, theta, fc='black', ec='black'))
            # Add the path to the axis
            self.ax.add_patch(self.handle[1])
        else: self.handle[1].set_angle(np.rad2deg(theta))
