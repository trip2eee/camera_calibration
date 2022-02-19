"""
@fn     generate_test_data.py
@brief  This file generates calibration test data.
@author Jongmin Park
@date   February 02, 2022
"""

from homography import Homography
import numpy as np
import matplotlib.pyplot as plt

param = {
    'image_width':1280,
    'image_height':800,
    'fu':1200,
    'fv':1200,
    'cu':640,
    'cv':400,
    'dist':[-1, -0.5, -0.9, 0.0, 0.5, 1.0],

    'roll': -0.03,
    'pitch': 0.02,
    'yaw': 0.01,
    'tx': 1.0,
    'ty': 0.1,
    'tz': -1.0
}

target = {
    'x':1.5,
    'z':0.2,
    'y':-1,
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

    x = np.linspace(x0, x0, 1)
    y = np.linspace(y0, y0+(cols-1)*check_size, cols)
    z = np.linspace(z0, z0+(rows-1)*check_size, rows)

    points = np.meshgrid(x, y, z)

    target_x = points[0].flatten()
    target_y = points[1].flatten()
    target_z = points[2].flatten()

    return target_x, target_y, target_z

if __name__ == '__main__':
    
    homography = Homography(**param)
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



