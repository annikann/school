"""
Class for plotting Assignment 1 Sample UAV

"""

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


class UAV_animation():
    def __init__(self, state0, scale=0.25, multfigs=False):

        self.scale = scale
        self.flag_init = True
        self.fig = plt.figure(1)
        if multfigs == True:
            self.ax = self.fig.add_subplot(1, 2, 1, projection="3d")
        else:
            self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.set_xlim([-10,10])
        self.ax.set_ylim([-10,10])
        self.ax.set_zlim([-10,10])
        self.ax.set_title('3D Animation')
        self.ax.set_xlabel('East(m)')
        self.ax.set_ylabel('North(m)')
        self.ax.set_zlabel('Height(m)')
        
        #self.update(state0)
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
                            [fuse_l2, -fuse_w/2, fuse_h/2],
                            [fuse_l2, -fuse_w/2, -fuse_h/2],
                            [fuse_l2, fuse_w/2, -fuse_h/2],
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
                            [-fuse_l3, 0, vtail_h]])
    

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
                            [vr[5], vr[14], vr[15]]])*[1, 1, -1]

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
            poly = Poly3DCollection(verts, facecolors=['g', 'r', 'r', 'r', 'r','r'], alpha=.6)
            self.UAV = self.ax.add_collection3d(poly)#
            plt.pause(0.01)
        else:
            self.UAV.set_verts(verts)
            plt.pause(0.01)

    def transformation_matrix(self):
        x = self.x
        y = self.y
        z = self.z
        roll = -self.roll
        pitch = -self.pitch
        yaw = self.yaw
        return np.array(
            [[cos(yaw) * cos(pitch), -sin(yaw) * cos(roll) + cos(yaw) * sin(pitch) * sin(roll), sin(yaw) * sin(roll) + cos(yaw) * sin(pitch) * cos(roll), x],
             [sin(yaw) * cos(pitch), cos(yaw) * cos(roll) + sin(yaw) * sin(pitch)
              * sin(roll), -cos(yaw) * sin(roll) + sin(yaw) * sin(pitch) * cos(roll), y],
             [-sin(pitch), cos(pitch) * sin(roll), cos(pitch) * cos(yaw), z]
             ])

    def plot(self):  # pragma: no cover
        T = self.transformation_matrix()

        p1_t = np.matmul(T, self.p1)
        p2_t = np.matmul(T, self.p2)
        p3_t = np.matmul(T, self.p3)
        p4_t = np.matmul(T, self.p4)

        #plt.cla() # use handle 
        if self.flag_init is True:
            body, =self.ax.plot([p1_t[0], p2_t[0], p3_t[0], p4_t[0], p1_t[0], p3_t[0], p4_t[0], p2_t[0]],
                        [p1_t[1], p2_t[1], p3_t[1], p4_t[1], p1_t[1], p3_t[1], p4_t[1], p2_t[1]],
                        [p1_t[2], p2_t[2], p3_t[2], p4_t[2], p1_t[2], p3_t[2], p4_t[2], p2_t[2]], 'k-') # rotor
            self.handle.append(body)

            traj, =self.ax.plot(self.x_data, self.y_data, self.z_data, 'b:')# trajectory
            self.handle.append(traj)

            plt.xlim(-2, 2)
            plt.ylim(-2, 2)
            self.ax.set_zlim(0, 4)
            plt.xlabel('North')
            plt.ylabel('East')
            self.flag_init = False 
            plt.pause(0.001) # can be put in the main file
        else:
            self.handle[0].set_data([p1_t[0], p2_t[0], p3_t[0], p4_t[0], p1_t[0], p3_t[0], p4_t[0], p2_t[0]],
                        [p1_t[1], p2_t[1], p3_t[1], p4_t[1], p1_t[1], p3_t[1], p4_t[1], p2_t[1]])
            self.handle[0].set_3d_properties([p1_t[2], p2_t[2], p3_t[2], p4_t[2],p1_t[2], p3_t[2], p4_t[2], p2_t[2]])


            self.handle[1].set_data(self.x_data, self.y_data)
            self.handle[1].set_3d_properties(self.z_data)
            print(self.handle)
            plt.pause(0.001)
