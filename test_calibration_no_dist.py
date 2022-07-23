"""
@fn     test_calibration_no_dist.py
@brief  This file test camera calibration with data without distortion.
@author Jongmin Park
@date   July 23, 2022
"""
from camera_geometry import CameraGeometry
from generate_test_data import generate_target_points
from calibrate_camera import CalibrateCamera

import matplotlib.pyplot as plt
import numpy as np

param1 = {
    'image_width':1280,
    'image_height':800,
    'fu':1000,
    'fv':1000,
    'cu':640,
    'cv':400,
    'dist':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    'roll': 0.0,
    'pitch': 0.2,
    'yaw': 0.0,
    'tx': -0.1,
    'ty': 0.0,
    'tz': 3.0
}

param2 = {
    'image_width':1280,
    'image_height':800,
    'fu':1000,
    'fv':1000,
    'cu':640,
    'cv':400,
    'dist':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    'roll': 0.0,
    'pitch': 0.0,
    'yaw': 0.05,
    'tx': 0.0,
    'ty': 0.1,
    'tz': 3.0
}

param3 = {
    'image_width':1280,
    'image_height':800,
    'fu':1000,
    'fv':1000,
    'cu':640,
    'cv':400,
    'dist':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    'roll': 0.02,
    'pitch': -0.1,
    'yaw': 0.1,
    'tx': 0.1,
    'ty': 0.0,
    'tz': 3.2
}

param4 = {
    'image_width':1280,
    'image_height':800,
    'fu':1000,
    'fv':1000,
    'cu':640,
    'cv':400,
    'dist':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
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


if __name__ == '__main__':
    
    target_x, target_y, target_z = generate_target_points(**target)
    pt_world = np.vstack((target_x, target_y, target_z)).T

    homography = CameraGeometry(**param1)
    target_u, target_v = homography.xyz_to_uv(target_x, target_y, target_z)    
    pt_image1 = np.vstack((target_u, target_v)).T

    homography = CameraGeometry(**param2)
    target_u, target_v = homography.xyz_to_uv(target_x, target_y, target_z)
    pt_image2 = np.vstack((target_u, target_v)).T

    homography = CameraGeometry(**param3)
    target_u, target_v = homography.xyz_to_uv(target_x, target_y, target_z)
    pt_image3 = np.vstack((target_u, target_v)).T

    homography = CameraGeometry(**param4)
    target_u, target_v = homography.xyz_to_uv(target_x, target_y, target_z)
    pt_image4 = np.vstack((target_u, target_v)).T

    list_pt_image = [pt_image1, pt_image2, pt_image3, pt_image4]

    calib = CalibrateCamera()
    calib.calibrate(pt_world, list_pt_image)


    # verify Homography

    for idx_image, pt_image in enumerate(list_pt_image):
        x0 = calib.target_origin[0]
        y0 = calib.target_origin[1]
        H = calib.homography[idx_image]

        for i in range(len(pt_world)):
            x = pt_world[i,0] - x0
            y = pt_world[i,1] - y0
            xy1 = np.matrix([[x], [y], [1]])
            p = np.matmul(H, xy1)
            u = p[0,0]/p[2,0]
            v = p[1,0]/p[2,0]

            err_u = pt_image[i,0] - u
            err_v = pt_image[i,1] - v
            sqr_err = err_u**2 + err_v**2

            if sqr_err > 1e-20:
                print('Inaccurate Homography')
                break

    for idx_image, pt_image in enumerate(list_pt_image):
        plt.figure('image {}'.format(idx_image))
        plt.plot([0, param1['image_width']], [0, 0])
        plt.plot([0, param1['image_width']], [param1['image_height'], param1['image_height']])
        plt.plot([0, 0], [0, param1['image_height']])
        plt.plot([param1['image_width'], param1['image_width']], [0, param1['image_height']])

        plt.scatter(pt_image[:,0], pt_image[:,1])
        plt.gca().invert_yaxis()
        plt.axis('equal')
    plt.show()



