"""
@fn     visualizer.py
@brief  This module visualizes calibration results.
        repository: https://github.com/trip2eee/camera_calibration
@author Jongmin Park
@date   July 25, 2022
"""

import numpy as np
from camera_geometry import CameraGeometry
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Visualizer:
    def __init__(self):
        pass

    def compute_RT(self, extrinsic):
        roll  = extrinsic['roll']
        pitch = extrinsic['pitch']        
        yaw   = extrinsic['yaw']

        tx = extrinsic['tx']
        ty = extrinsic['ty']
        tz = extrinsic['tz']

        cosr = np.cos(roll)
        sinr = np.sin(roll)
        Rx = np.array([[1, 0,     0,  ],
                       [0, cosr, -sinr],
                       [0, sinr,  cosr]])

        # pitch
        cosp = np.cos(pitch)
        sinp = np.sin(pitch)
        Ry = np.array([[cosp,  0, sinp],
                       [0,     1,    0],
                       [-sinp, 0, cosp]])

        # yaw
        cosy = np.cos(yaw)
        siny = np.sin(yaw)
        Rz = np.array([[cosy, -siny, 0],
                       [siny,  cosy, 0],
                       [0,        0, 1]])
        R  = np.matmul(Rz, np.matmul(Ry, Rx))

        # translation
        T = np.array([[tx],
                      [ty],
                      [tz]])        
        T = np.matmul(R, T)

        return R, T

    def draw_axis(self, ax:Axes3D):
        ax.plot([0, 1],[0, 0],[0, 0], c='r')
        ax.plot([0, 0],[0, 1],[0, 0], c='g')
        ax.plot([0, 0],[0, 0],[0, 1], c='b')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

    def draw_target_plane(self, ax:Axes3D, pt1, pt2, R, T, c):
        # p_c = [R RT]p1
        # R^-1*p_c = p + T
        # p = R^-1*p_c - T

        x1, y1 = pt1
        x2, y2 = pt2

        # (x1, y1)      (x2, y1)
        #     +------------+
        #     |            |
        #     +------------+
        # (x1, y2)      (x2, y2)

        # left top
        p_lt = np.array([[x1], [y1], [0]])
        # right top
        p_rt = np.array([[x2], [y1], [0]])
        # left bottom
        p_lb = np.array([[x1], [y2], [0]])
        # right bottom
        p_rb = np.array([[x2], [y2], [0]])

        p_lt = np.matmul(R.T, p_lt) - T
        p_rt = np.matmul(R.T, p_rt) - T
        p_lb = np.matmul(R.T, p_lb) - T
        p_rb = np.matmul(R.T, p_rb) - T

        ax.plot([p_lt[0,0], p_rt[0,0],  p_rb[0,0],  p_lb[0,0],  p_lt[0,0]], 
                [p_lt[1,0], p_rt[1,0],  p_rb[1,0],  p_lb[1,0],  p_lt[1,0]],
                [p_lt[2,0], p_rt[2,0],  p_rb[2,0],  p_lb[2,0],  p_lt[2,0]], c=c)

        ax.text(p_lb[0,0], p_lb[1,0], p_lb[2,0], 'target', color=c)

    def draw_camera(self, ax, R, T, c, cam_id):
        # (x1, y1)      (x2, y1)
        #     +------------+
        #     |            |
        #     +------------+
        # (x1, y2)      (x2, y2)

        # p_c = [R RT]p1
        # R^-1*p_c = p + T
        # p = R^-1*p_c - T

        s = 0.25
        x1, y1 = -s, -s
        x2, y2 = s, s
        z = s*2

        # left top
        p_lt = np.array([[x1], [y1], [z]])
        # right top
        p_rt = np.array([[x2], [y1], [z]])
        # left bottom
        p_lb = np.array([[x1], [y2], [z]])
        # right bottom
        p_rb = np.array([[x2], [y2], [z]])
        # center
        p_ct = np.array([[0], [0], [0]])

        p_lt = np.matmul(R.T, p_lt) - T
        p_rt = np.matmul(R.T, p_rt) - T
        p_lb = np.matmul(R.T, p_lb) - T
        p_rb = np.matmul(R.T, p_rb) - T
        p_ct = np.matmul(R.T, p_ct) - T

        label = 'cam{}'.format(cam_id)
        # (lt - rt - rb - lb - lt) - (ct - rt - rb - ct - lb)
        ax.plot([p_lt[0,0], p_rt[0,0], p_rb[0,0], p_lb[0,0], p_lt[0,0], p_ct[0,0], p_rt[0,0], p_rb[0,0], p_ct[0,0], p_lb[0,0]], 
                [p_lt[1,0], p_rt[1,0], p_rb[1,0], p_lb[1,0], p_lt[1,0], p_ct[1,0], p_rt[1,0], p_rb[1,0], p_ct[1,0], p_lb[1,0]],
                [p_lt[2,0], p_rt[2,0], p_rb[2,0], p_lb[2,0], p_lt[2,0], p_ct[2,0], p_rt[2,0], p_rb[2,0], p_ct[2,0], p_lb[2,0]], c=c, label=label)

        ax.text(p_ct[0,0], p_ct[1,0], p_ct[2,0] - 0.1, label, color=c)

    def visualize_cameras(self, extrinsics):
        """This method visualizes camera(positions and angles) in the world coordinate.
        """
        fig = plt.figure('Cameras')        
        ax = fig.add_subplot(projection='3d')

        self.draw_axis(ax)

        range = 4
        ax.set_xlim(-range*0.5, range*0.5)
        ax.set_ylim(-range*0.5, range*0.5)
        ax.set_zlim(-range, 0)

        # draw target
        tw = 1
        th = 1

        pt1 = [-tw, -th]
        pt2 = [+tw, +th]
        R = np.eye(3)
        T = np.zeros((3, 1))

        self.draw_target_plane(ax, pt1, pt2, R, T, 'k')

        color_list = ['r', 'g', 'b', 'c', 'm', 'y']
        for idx_camera, extrinsic in enumerate(extrinsics):
            R, T = self.compute_RT(extrinsic)
            self.draw_camera(ax, R, T, color_list[idx_camera % 6], idx_camera)
        plt.legend()
        plt.show()

    def visualize_reprojection(self, intrinsic, extrinsics, pt_world, pt_images):

        for idx_image, extrinsic in enumerate(extrinsics):
            param = extrinsic
            param['intrinsic'] = intrinsic
            cam_geometry = CameraGeometry(**param)

            u = pt_images[idx_image][:,0]
            v = pt_images[idx_image][:,1]
            u_reproj, v_reproj = cam_geometry.xyz_to_uv(pt_world[:,0], pt_world[:,1], pt_world[:,2])

            plt.figure('image{}'.format(idx_image))
            plt.scatter(u, v, c='g', label='target')
            plt.scatter(u_reproj, v_reproj, c='b', label='reprojection')
            plt.axis('equal')
            plt.legend()
        
        plt.show()


