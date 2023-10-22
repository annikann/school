# """
# Class for plotting Assignment 1 Sample UAV

# """

import sys
import os
sys.path.append('/Users/annikacarlson/Documents/FlightMech/flight_sim')
from math import cos, sin
import numpy as np
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from mpl_toolkits.mplot3d import Axes3D
from library.rotations import Quaternion2Euler, Quaternion2Rotation, Euler2Rotation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

class MAV_animation:
    def __init__(self, limits=10, scale=0.25, multfigs=False):

        self.scale = scale
        self.flag_init = True
        self.fig = plt.figure(1)
        if multfigs == True:
            self.ax = self.fig.add_subplot(1, 3, 2, projection="3d")
        else:
            self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.set_xlim([-limits,limits])
        self.ax.set_ylim([-limits,limits])
        self.ax.set_zlim([-limits,limits])
        self.lim = limits
        self.ax.set_title('3D Animation')
        self.ax.set_xlabel('East(m)')
        self.ax.set_ylabel('North(m)')
        self.ax.set_zlabel('Height(m)')

    def UAV_vertices(self,pn,pe,pd,phi,theta,psi):
        
        fuse_l1 = 2.8
        fuse_l2 = 0.1
        fuse_l3 = 5
        fuse_h = 1
        fuse_w = 1
        wing_l = 1.5
        wing_w = 7
        vtail_h = 1.5
        vtail_l = 1
        htail_l = 0.75
        htail_w = 3.3

        UAV_verts = np.array([[fuse_l1, 0, 0],
                            [fuse_l2, fuse_w/2, fuse_h/2],
                            [fuse_l2, fuse_w/2, -fuse_h/2],
                            [fuse_l2, -fuse_w/2, fuse_h/2],
                            [fuse_l2, -fuse_w/2, -fuse_h/2],
                            [-fuse_l3, 0, 0],
                            [0, wing_w/2, 0],
                            [-wing_l, wing_w/2, 0],
                            [-wing_l, -wing_w/2, 0],
                            [0, -wing_w/2, 0],
                            [-(fuse_l3 - htail_l), htail_w/2, 0],
                            [-fuse_l3, htail_w/2, 0],
                            [-fuse_l3, -htail_w/2, 0],
                            [-(fuse_l3 - htail_l), -htail_w/2, 0],
                            [-(fuse_l3 - vtail_l), 0, 0],
                            [-fuse_l3, 0, -vtail_h]])
    

        pos_ned=np.array([pn, pe, pd])

        # create m by n copies of pos_ned and used for translation
        ned_rep= np.tile(pos_ned.copy(), (np.shape(UAV_verts)[0], 1))

        R = Euler2Rotation(phi,theta,psi)
        vr = np.matmul(R,UAV_verts.T).T
        vr = vr + ned_rep

        # rotate for plotting north=y east=x h=-z
        R_plot = np.array([[0, 1, 0],
                           [1, 0, 0],
                           [0, 0, -1]])
        
        vr = np.matmul(R_plot,vr.T).T

        UAV_faces = np.array([[vr[0], vr[1], vr[2]],
                            [vr[0], vr[1], vr[4]],
                            [vr[0], vr[2], vr[3]],
                            [vr[1], vr[3], vr[4]],
                            [vr[2], vr[3], vr[5]],
                            [vr[1], vr[2], vr[5]],
                            [vr[1], vr[4], vr[5]],
                            [vr[3], vr[4], vr[5]],
                            [vr[6], vr[7], vr[8]],
                            [vr[6], vr[8], vr[9]],
                            [vr[10], vr[11], vr[12]],
                            [vr[10], vr[12], vr[13]],
                            [vr[5], vr[14], vr[15]]])*[1, 1, 1]

        return(UAV_faces)

    def update(self, pn, pe, pd, phi, theta, psi):
        # draw plot elements
        self.draw_UAV(pn,pe,pd,phi,theta,psi)
        # Set initialization flag to False after first call
        if self.flag_init == True:
            self.flag_init = False

    def draw_UAV(self, pn, pe, pd, phi, theta, psi):
        verts = self.UAV_vertices(pn,pe,pd,phi,theta,psi)
        if self.flag_init is True:
            poly = Poly3DCollection(verts, facecolors=["b"], alpha=.6)
            self.UAV = self.ax.add_collection3d(poly)#
            self.ax.set_xlim([pe-self.lim, pe+self.lim])
            self.ax.set_ylim([pn-self.lim, pn+self.lim])
            self.ax.set_zlim([-pd-self.lim, -pd+self.lim])
            plt.pause(0.01)
        else:
            self.UAV.set_verts(verts)
            self.ax.set_xlim([pe-self.lim, pe+self.lim])
            self.ax.set_ylim([pn-self.lim, pn+self.lim])
            self.ax.set_zlim([-pd-self.lim, -pd+self.lim])
            plt.pause(0.01)
