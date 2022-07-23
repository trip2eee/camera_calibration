"""
@fn     calibrate_camera.py
@brief  This file calibrates camera
@author Jongmin Park
@date   July 23, 2022
"""

import numpy as np


class CalibrateCamera:
    def __init__(self):
        self.target_origin = None
        self.num_images = 0
        self.homography = []

        # intrinsic camera parameters
        self.cu = 0
        self.cv = 0
        self.fu = 0
        self.fv = 0

    def calibrate(self, pt_world, list_pt_image):
        self.num_images = len(list_pt_image)

        self.target_origin = pt_world[0,:]
        
        for idx_image in range(len(list_pt_image)):
            self.compute_homography(pt_world, list_pt_image[idx_image])
        
        self.compute_initial_value()

    def compute_initial_value(self):
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
            

            Ai = np.matrix([[h_11*h_12, h_11*h_32 + h_12*h_31, h_21*h_22, h_21*h_32 + h_22*h_31, h_31*h_32], 
                            [h_11**2 - h_12**2, 2*h_11*h_31 - 2*h_12*h_32, h_21**2 - h_22**2, 2*h_21*h_31 - 2*h_22*h_32, h_31**2 - h_32**2]])
        
            A[idx_image*2:idx_image*2+2, :] = Ai

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
        
        self.cu = cu
        self.cv = cv
        self.fu = fu
        self.fv = fv

    def compute_homography(self, pt_world, pt_image):
        """This method computes homography matrix for the given pair of world coordinates and image coordinates
        """
        assert(len(pt_world) == len(pt_image))

        A = np.zeros([8, 8])
        b = np.zeros([8, 1])

        x0 = self.target_origin[0]
        y0 = self.target_origin[1]
        
        for i in range(len(pt_world)):
            x = pt_world[i,0] - x0
            y = pt_world[i,1] - y0
            # z = 0

            u = pt_image[i,0]
            v = pt_image[i,1]

            Ai = np.matrix([[-x, -y, -1,  0,  0,  0, u*x, u*y],
                            [ 0,  0,  0, -x, -y, -1, v*x, v*y]])
            bi = np.matrix([[-u], 
                            [-v]])
            
            A += Ai.T@Ai
            b += Ai.T@bi


        h = np.linalg.solve(A, b)

        H = np.matrix([[h[0,0], h[1,0], h[2,0]],
                       [h[3,0], h[4,0], h[5,0]],
                       [h[6,0], h[7,0], 1.0]])
        
        self.homography.append(H)



        
