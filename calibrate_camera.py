"""
@fn     calibrate_camera.py
@brief  This file calibrates camera
@author Jongmin Park
@date   July 23, 2022
"""

from re import S
import numpy as np


class CalibrateCamera:
    def __init__(self):
        self.target_origin = None
        self.num_images = 0
        
        # intrinsic camera parameters        
        self.intrinsic = {
            'image_width':0,
            'image_height':0,
            'fu':0,
            'fv':0,
            'cu':0,
            'cv':0,
            'dist':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        }

        self.homography = []
        self.extrinsics = []

    def calibrate(self, pt_world, list_pt_image):
        self.num_images = len(list_pt_image)
        

        self.target_origin = pt_world[0,:]
        
        for idx_image in range(len(list_pt_image)):
            self.compute_homography(pt_world, list_pt_image[idx_image])
        
        self.compute_initial_values()

    def compute_initial_values(self):
        """This method computes initial parameters
        """
        A = np.zeros((2 * self.num_images, 5))

        for idx_image, H in enumerate(self.homography):
            h_11 = H[0,0]
            h_21 = H[1,0]
            h_31 = H[2,0]

            h_12 = H[0,1]
            h_22 = H[1,1]
            h_32 = H[2,1]

            A[idx_image*2:idx_image*2+2, :] = np.array([[h_11*h_12, h_11*h_32 + h_12*h_31, h_21*h_22, h_21*h_32 + h_22*h_31, h_31*h_32], 
                                                        [h_11**2 - h_12**2, 2*h_11*h_31 - 2*h_12*h_32, h_21**2 - h_22**2, 2*h_21*h_31 - 2*h_22*h_32, h_31**2 - h_32**2]])

        u, d, v = np.linalg.svd(A)
        v = v.T
        b = v[:,-1]

        # determine sign
        if b[0] < 0.0:
            b *= -1.0

        b11 = b[0]
        b13 = b[1]
        b22 = b[2]
        b23 = b[3]
        b33 = b[4]

        cu = -b13 / b11
        cv = -b23 / b22

        # scale: 1/lambda in README.md
        scale = b33 + b13*cu + b23*cv

        fu = 1 / np.sqrt(b11 / scale)
        fv = 1 / np.sqrt(b22 / scale)
        
        self.intrinsic['fu'] = fu
        self.intrinsic['fv'] = fv
        self.intrinsic['cu'] = cu
        self.intrinsic['cv'] = cv
                
        K = np.array([[fu,  0, cu],
                      [ 0, fv, cv],
                      [ 0,  0,  1]])
        invK = np.linalg.inv(K)

        # compute extrinsic camera parameters
        for idx_image, H in enumerate(self.homography):
            sRT = np.matmul(invK, H)

            sr1 = sRT[:,0]
            sr2 = sRT[:,1]
            srt = sRT[:,2]

            # scale factor
            scale = np.sqrt(np.matmul(sr1.T, sr1))

            r1 = sr1 / scale
            r2 = sr2 / scale
            r3 = np.cross(r1, r2)

            rt = srt / scale

            r1 = r1.reshape([3,1])
            r2 = r2.reshape([3,1])
            r3 = r3.reshape([3,1])
            rt = rt.reshape([3,1])

            R = np.hstack([r1, r2, r3])
            t = np.matmul(R.T, rt)
            
            pitch = -np.arcsin(R[2,0])
            yaw = np.arctan2(R[1,0], R[0,0])
            roll = np.arctan2(R[2,1], R[2,2])

            tx = t[0,0]
            ty = t[1,0]
            tz = t[2,0]

            extrinsic = {
                'roll': roll,
                'pitch': pitch,
                'yaw': yaw,
                'tx': tx,
                'ty': ty,
                'tz': tz
            }
            self.extrinsics.append(extrinsic)

    def compute_homography(self, pt_world, pt_image):
        """This method computes homography matrix for the given pair of world coordinates and image coordinates
        """
        assert(len(pt_world) == len(pt_image))

        A = np.zeros([8, 8])
        b = np.zeros([8, 1])

        for i in range(len(pt_world)):
            x = pt_world[i,0]
            y = pt_world[i,1]
            # z = 0

            u = pt_image[i,0]
            v = pt_image[i,1]

            Ai = np.array([[-x, -y, -1,  0,  0,  0, u*x, u*y],
                           [ 0,  0,  0, -x, -y, -1, v*x, v*y]])
            bi = np.array([[-u], 
                           [-v]])
            
            A += np.matmul(Ai.T, Ai)
            b += np.matmul(Ai.T, bi)


        h = np.linalg.solve(A, b)

        H = np.array([[h[0,0], h[1,0], h[2,0]],
                      [h[3,0], h[4,0], h[5,0]],
                      [h[6,0], h[7,0], 1.0]])
        
        self.homography.append(H)



        
