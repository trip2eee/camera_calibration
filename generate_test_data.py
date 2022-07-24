"""
@fn     generate_test_data.py
@brief  This file generates calibration test data.
@author Jongmin Park
@date   February 02, 2022
"""

from camera_geometry import CameraGeometry
import numpy as np
import matplotlib.pyplot as plt

param = {
    'image_width':1280,
    'image_height':800,
    'fu':1000,
    'fv':1000,

    'cu':640,
    'cv':400,
    # 'dist':[-0.3, -0.03, -0.003, 0.0, 0.0, 0.0],
    'dist':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],

    'roll': 0.0,
    'pitch': 0.0,
    'yaw': 0.0,
    'tx': 0.0,
    'ty': 0.0,
    'tz': 0.0
}

target = {
    'x':-1.0,
    'y':-1.0,
    'z':3,
    
    'check_size':0.2,
    'rows':10,
    'cols':10
}

def generate_target_points(**target):
    rows = target['rows']
    cols = target['cols']
    check_size = target['check_size']
    x0 = target['x']
    y0 = target['y']
    z0 = target['z']

    x = np.linspace(x0, x0+(cols-1)*check_size, cols)
    y = np.linspace(y0, y0+(rows-1)*check_size, rows)
    z = np.linspace(z0, z0,1)

    points = np.meshgrid(x, y, z)

    target_x = points[0].flatten()
    target_y = points[1].flatten()
    target_z = points[2].flatten()

    pt_world = np.vstack((target_x, target_y, target_z)).T

    return pt_world

if __name__ == '__main__':
    
    homography = CameraGeometry(**param)
    target_x, target_y, target_z = generate_target_points(**target)
    target_u, target_v = homography.xyz_to_uv(target_x, target_y, target_z)
    
    plt.figure('image')
    plt.plot([0, param['image_width']], [0, 0])
    plt.plot([0, param['image_width']], [param['image_height'], param['image_height']])
    plt.plot([0, 0], [0, param['image_height']])
    plt.plot([param['image_width'], param['image_width']], [0, param['image_height']])

    plt.scatter(target_u, target_v)
    plt.gca().invert_yaxis()
    plt.axis('equal')
    plt.show()



