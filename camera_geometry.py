"""
@fn     camera_geometry.py
@brief  CameraGeometry class.
        repository: https://github.com/trip2eee/camera_calibration
@author Jongmin Park
@date   February 02, 2022
"""

import numpy as np

class CameraGeometry:
    """ Homography class
    """
    def __init__(self, **params):
        """This constructor reads parameters and calls computeH() for computing homography matrix from the given parameters.
        """
        intrinsic = params['intrinsic']

        self.image_width = intrinsic['image_width']
        self.image_height = intrinsic['image_height']
        self.fu = intrinsic['fu']
        self.fv = intrinsic['fv']
        self.cu = intrinsic['cu']
        self.cv = intrinsic['cv']
        self.dist = intrinsic['dist']  # distortion parameters

        self.roll = params['roll']
        self.pitch = params['pitch']
        self.yaw = params['yaw']
        self.tx = params['tx']
        self.ty = params['ty']
        self.tz = params['tz']
        
        self.computeH()

    def computeH(self):
        """This method computes homography matrix from the given parameters.
        """
        cosr = np.cos(self.roll)
        sinr = np.sin(self.roll)
        Rx = np.array([[1, 0,     0,  ],
                       [0, cosr, -sinr],
                       [0, sinr,  cosr]])

        # pitch
        cosp = np.cos(self.pitch)
        sinp = np.sin(self.pitch)
        Ry = np.array([[cosp,  0, sinp],
                       [0,     1,    0],
                       [-sinp, 0, cosp]])

        # yaw
        cosy = np.cos(self.yaw)
        siny = np.sin(self.yaw)
        Rz = np.array([[cosy, -siny, 0],
                       [siny,  cosy, 0],
                       [0,        0, 1]])

        # translation
        T = np.array([[1, 0, 0, self.tx],
                      [0, 1, 0, self.ty],
                      [0, 0, 1, self.tz]])

        # intrinsic camera parameter matrix
        fu = self.fu
        fv = self.fv
        cu = self.cu
        cv = self.cv
        K = np.array([[fu, 0,  cu],
                      [0,  fv, cv],
                      [0,  0,  1,]])
        
        R  = np.matmul(Rz, np.matmul(Ry, Rx))

        self.R = R
        self.T = T
        
        RT = np.matmul(R, T)
        H  = np.matmul(K, RT)
        
        self.H = H
        self.H3x3 = np.array([[H[0,0], H[0,1], H[0,3]],
                              [H[1,0], H[1,1], H[1,3]],
                              [H[2,0], H[2,1], H[2,3]]])
        self.invH3x3 = np.linalg.inv(self.H3x3)

        
    def xyz_to_uv(self, x:np.ndarray, y:np.ndarray, z:np.ndarray):

        len_x = len(x)
        ones = np.ones(len_x)
        xyz = np.vstack((x, y, z, ones))
        
        # project xyz into image (undistorted).
        p = np.matmul(self.H, xyz)
        
        u = p[0,:] / p[2,:]
        v = p[1,:] / p[2,:]

        u = np.squeeze(np.asarray(u))
        v = np.squeeze(np.asarray(v))
        
        # normalize points.
        u0 = (u - self.cu) / self.fu
        v0 = (v - self.cv) / self.fv

        # distort image points.
        r2 = np.power(u0, 2) + np.power(v0, 2)
        r4 = np.power(r2, 2)
        r6 = np.power(r2, 3)

        k1 = self.dist[0]
        k2 = self.dist[1]
        k3 = self.dist[2]
        k4 = self.dist[3]
        k5 = self.dist[4]
        k6 = self.dist[5]

        ratios = (1.0 + k1*r2 + k2*r4 + k3*r6) / (1.0 + k4*r2 + k5*r4 + k6*r6)

        u = self.fu * np.multiply(ratios, u0) + self.cu
        v = self.fv * np.multiply(ratios, v0) + self.cv

        return u, v

