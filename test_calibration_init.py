"""
@fn     test_calibration_init.py
@brief  This file tests computation of initial values with simulated data.
        - No distortion is considered.
@author Jongmin Park
@date   July 23, 2022
"""
from camera_geometry import CameraGeometry
from generate_test_data import generate_target_points
from calibrate_camera import CalibrateCamera

import matplotlib.pyplot as plt
import numpy as np
import unittest

intrinsic = {
    'image_width':1280,
    'image_height':800,
    'fu':1000,
    'fv':1000,
    'cu':640,
    'cv':400,
    'dist':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
}

param1 = {
    'intrinsic':intrinsic,
    'roll': 0.0,
    'pitch': 0.2,
    'yaw': 0.0,
    'tx': -0.1,
    'ty': 0.0,
    'tz': 3.0
}

param2 = {
    'intrinsic':intrinsic,
    'roll': 0.0,
    'pitch': 0.0,
    'yaw': 0.05,
    'tx': 0.0,
    'ty': 0.1,
    'tz': 3.0
}

param3 = {
    'intrinsic':intrinsic,
    'roll': 0.02,
    'pitch': -0.1,
    'yaw': 0.1,
    'tx': 0.1,
    'ty': 0.0,
    'tz': 3.2
}

param4 = {
    'intrinsic':intrinsic,
    'roll': 0.02,
    'pitch': +0.1,
    'yaw': -0.1,
    'tx': 0.1,
    'ty': 0.4,
    'tz': 3.2
}

target = {
    'x':-1.0,
    'y':-1.0,
    'z':0,
    'check_size':0.2,
    'rows':10,
    'cols':10
}

params = [
    param1,
    param2,
    param3,
    param4
]

class TestCalibration_NoDist(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_calibrate(self):
        pt_world = generate_target_points(**target)        

        homography = CameraGeometry(**param1)
        target_u, target_v = homography.xyz_to_uv(pt_world[:,0], pt_world[:,1], pt_world[:,2])
        pt_image1 = np.vstack((target_u, target_v)).T

        homography = CameraGeometry(**param2)
        target_u, target_v = homography.xyz_to_uv(pt_world[:,0], pt_world[:,1], pt_world[:,2])
        pt_image2 = np.vstack((target_u, target_v)).T

        homography = CameraGeometry(**param3)
        target_u, target_v = homography.xyz_to_uv(pt_world[:,0], pt_world[:,1], pt_world[:,2])
        pt_image3 = np.vstack((target_u, target_v)).T

        homography = CameraGeometry(**param4)
        target_u, target_v = homography.xyz_to_uv(pt_world[:,0], pt_world[:,1], pt_world[:,2])
        pt_image4 = np.vstack((target_u, target_v)).T

        list_pt_image = [pt_image1, pt_image2, pt_image3, pt_image4]

        calib = CalibrateCamera()
        calib.calibrate(pt_world, list_pt_image)

        # test for homography
        for idx_image, pt_image in enumerate(list_pt_image):            
            H = calib.homography[idx_image]

            for i in range(len(pt_world)):
                x = pt_world[i,0]
                y = pt_world[i,1]
                xy1 = np.array([[x], [y], [1]])
                p = np.matmul(H, xy1)
                u = p[0,0]/p[2,0]
                v = p[1,0]/p[2,0]

                err_u = pt_image[i,0] - u
                err_v = pt_image[i,1] - v
                sqr_err = err_u**2 + err_v**2

                self.assertLess(sqr_err, 1e-20, 'Inaccurate Homography')

        # test for intrinsic
        self.assertAlmostEqual(calib.intrinsic['fu'], intrinsic['fu'], 6)
        self.assertAlmostEqual(calib.intrinsic['fv'], intrinsic['fv'], 6)
        self.assertAlmostEqual(calib.intrinsic['cu'], intrinsic['cu'], 6)
        self.assertAlmostEqual(calib.intrinsic['cv'], intrinsic['cv'], 6)
        
        # test for extrinsic
        for idx_image in range(calib.num_images):
            self.assertAlmostEqual(calib.extrinsics[idx_image]['roll'],  params[idx_image]['roll'], 6)
            self.assertAlmostEqual(calib.extrinsics[idx_image]['pitch'], params[idx_image]['pitch'], 6)
            self.assertAlmostEqual(calib.extrinsics[idx_image]['yaw'],   params[idx_image]['yaw'], 6)

            self.assertAlmostEqual(calib.extrinsics[idx_image]['tx'], params[idx_image]['tx'], 6)
            self.assertAlmostEqual(calib.extrinsics[idx_image]['ty'], params[idx_image]['ty'], 6)
            self.assertAlmostEqual(calib.extrinsics[idx_image]['tz'], params[idx_image]['tz'], 6)

        print(calib.intrinsic)

        image_width = intrinsic['image_width']
        image_height = intrinsic['image_height']
        for idx_image, pt_image in enumerate(list_pt_image):
            plt.figure('image {}'.format(idx_image))
            plt.plot([0, image_width], [0, 0])
            plt.plot([0, image_width], [image_height, image_height])
            plt.plot([0, 0], [0, image_height])
            plt.plot([image_width, image_width], [0, image_height])

            plt.scatter(pt_image[:,0], pt_image[:,1])
            plt.gca().invert_yaxis()
            plt.axis('equal')
        plt.show()

if __name__ == '__main__':
    
    unittest.main()
    

