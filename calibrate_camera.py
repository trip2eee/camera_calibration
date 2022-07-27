"""
@fn     calibrate_camera.py
@brief  This file calibrates camera
        repository: https://github.com/trip2eee/camera_calibration
@author Jongmin Park
@date   July 23, 2022
"""

import numpy as np
from quaternion import Quaternion
import time
from typing import List

MAX_ITERATIONS = 100                # maximum number of interation in refinement.
STOP_CONDITION_DELTA_PARAM = 1e-10  # minimum norm of delta param in refinement.
STOP_CONDITION_DAMPING = 1e5        # maximum damping factor in refinement (mu).

class CalibrateCamera:
    def __init__(self):        
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

        self.homography = []    # list of homography matrices for input images.
        self.extrinsics = []    # list of extrinsic camera parameters for input images.
        self.num_iterations = 20

    def calibrate(self, pt_world:np.ndarray, pt_images:List[np.ndarray], num_iterations=MAX_ITERATIONS):
        """This method is the main entry fuction.
           pt_world - World coordinates of corner points
           pt_images - List of image coordinates of corner points
           num_iterations - The number of iteration in refine_parameters().

           This method returns dictionary of intrinsic and list of dictionaries of extrinsics.
        """
        t_start = time.time()

        self.num_images= len(pt_images)
        self.num_iterations = num_iterations
        
        for idx_image in range(len(pt_images)):
            H = self.compute_homography(pt_world, pt_images[idx_image])
            self.homography.append(H)
        
        self.compute_initial_values()
        self.refine_parameters(pt_world, pt_images)

        t_end = time.time()
        print('Finished in {} seconds'.format(t_end - t_start))

        return self.intrinsic, self.extrinsics
    
    def print_parameters(self):
        """This method prints camera parameters
        """

        print('intrinsic parameters')
        print('fu, fv: {:.6f}, {:.6f}'.format(self.intrinsic['fu'], self.intrinsic['fv']))
        print('cu, cv: {:.6f}, {:.6f}'.format(self.intrinsic['cu'], self.intrinsic['cv']))
        print('distortion: ',end='')
        print(self.intrinsic['dist'])
        
        print('extrinsinc parameters (x, y, z)')
        for idx_image, extrinsic in enumerate(self.extrinsics):
            print('Angles {}: {:.6f}, {:.6f}, {:.6f}'.format(idx_image, extrinsic['roll'], extrinsic['pitch'], extrinsic['yaw']))
            print('Translation {}: {:.6f}, {:.6f}, {:.6f}'.format(idx_image, extrinsic['tx'], extrinsic['ty'], extrinsic['tz']))

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
           This method returns 3x3 homography matrix.
        """
        assert(len(pt_world) == len(pt_image))

        A = np.zeros([8, 8])
        b = np.zeros([8, 1])

        for i in range(len(pt_world)):
            x = pt_world[i,0]
            y = pt_world[i,1]

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
        
        return H

    def refine_parameters(self, pt_world, list_pt_image):
        """This method refines calibration parameters including distortion using Levenberg-Marquardt method.
        """
        num_intrinsics = 7
        num_extrinsics = 7
        num_images = len(list_pt_image)
        num_points = len(pt_world)

        num_rows_j = 3
        num_cols_j = num_intrinsics + (num_extrinsics * num_images)

        param = np.zeros([num_cols_j, 1])

        # initial values
        param[0,0] = self.intrinsic['fu']
        param[1,0] = self.intrinsic['fv']
        param[2,0] = self.intrinsic['cu']
        param[3,0] = self.intrinsic['cv']
        param[4,0] = self.intrinsic['dist'][0]
        param[5,0] = self.intrinsic['dist'][1]
        param[6,0] = self.intrinsic['dist'][2]

        for idx_image in range(num_images):
            roll  = self.extrinsics[idx_image]['roll']
            pitch = self.extrinsics[idx_image]['pitch']
            yaw   = self.extrinsics[idx_image]['yaw']
            tx    = self.extrinsics[idx_image]['tx']
            ty    = self.extrinsics[idx_image]['ty']
            tz    = self.extrinsics[idx_image]['tz']
            
            q_rot = Quaternion()
            q_rot.to_quaternion(roll, pitch, yaw)
            R = q_rot.to_rotation_matrix()
            t = np.array([[tx], [ty], [tz]])
            rt = np.matmul(R, t)

            param[num_intrinsics + (idx_image*num_extrinsics) + 0, 0] = q_rot.q[0]
            param[num_intrinsics + (idx_image*num_extrinsics) + 1, 0] = q_rot.q[1]
            param[num_intrinsics + (idx_image*num_extrinsics) + 2, 0] = q_rot.q[2]
            param[num_intrinsics + (idx_image*num_extrinsics) + 3, 0] = q_rot.q[3]
            param[num_intrinsics + (idx_image*num_extrinsics) + 4, 0] = rt[0,0]
            param[num_intrinsics + (idx_image*num_extrinsics) + 5, 0] = rt[1,0]
            param[num_intrinsics + (idx_image*num_extrinsics) + 6, 0] = rt[2,0]

        opt_ms_error = 0.0
        updateJ = True
        mu = 0.0001

        # compute initial error
        for idx_image in range(num_images):
            for idx_point in range(len(pt_world)):
                xyz = pt_world[idx_point]
                uv = list_pt_image[idx_image][idx_point]

                intrin = param[0:num_intrinsics, 0]
                extrin = param[num_intrinsics + (idx_image*num_extrinsics): num_intrinsics + ((idx_image+1)*num_extrinsics), 0]
                f = self.compute_f(intrin, extrin, uv, xyz)

                eu = f[0,0]
                ev = f[1,0]
                eq = f[2,0]
                opt_ms_error += eu**2 + ev**2 + eq**2
        opt_ms_error = opt_ms_error / (num_images * num_points)
        print('initial mean squared error: {}'.format(opt_ms_error))

        for iter in range(self.num_iterations):
            if updateJ:
                # JJ = J^T * J
                JJ = np.zeros([num_cols_j, num_cols_j])
                # Jf = J^T * f
                Jf = np.zeros([num_cols_j, 1])

                for idx_image in range(num_images):
                    for idx_point in range(len(pt_world)):
                        xyz = pt_world[idx_point]
                        uv = list_pt_image[idx_image][idx_point]

                        intrin = param[0:num_intrinsics, 0]
                        extrin = param[num_intrinsics + (idx_image*num_extrinsics): num_intrinsics + ((idx_image+1)*num_extrinsics), 0]
                        f = self.compute_f(intrin, extrin, uv, xyz)
                        J_intrin, J_extrin = self.compute_J(intrin, extrin, xyz)

                        J = np.zeros([num_rows_j, num_cols_j])
                        J[:,0:num_intrinsics] = J_intrin
                        J[:,num_intrinsics+idx_image*num_extrinsics:num_intrinsics+(idx_image+1)*num_extrinsics] = J_extrin

                        Jf += np.matmul(J.T, f)
                        JJ += np.matmul(J.T, J)
                updateJ = False

            # compute delta p
            for idx_d in range(num_cols_j):
                JJ[idx_d, idx_d] = (1.0 + mu) * JJ[idx_d, idx_d]
            
            delta_param = -np.linalg.solve(JJ, Jf)
            param_temp = param + delta_param

            # compute error of the temp. solution
            temp_ms_error = 0.0
            for idx_image in range(num_images):
                for idx_point in range(len(pt_world)):
                    xyz = pt_world[idx_point]
                    uv = list_pt_image[idx_image][idx_point]

                    intrin = param_temp[0:num_intrinsics, 0]
                    extrin = param_temp[num_intrinsics + (idx_image*num_extrinsics): num_intrinsics + ((idx_image+1)*num_extrinsics), 0]
                    f = self.compute_f(intrin, extrin, uv, xyz)
                    temp_ms_error += f[0,0]**2 + f[1,0]**2 + f[2,0]**2
            temp_ms_error = temp_ms_error / (num_images * num_points)
            
            if temp_ms_error < opt_ms_error:
                mu *= 0.1
                updateJ = True
                param = param_temp.copy()
                opt_ms_error = temp_ms_error
            else:
                mu *= 10

            print('iteration: {}, mean-squared projection error: {}'.format(iter, opt_ms_error))
            
            if mu >= STOP_CONDITION_DAMPING:
                break

            norm_delta_param = np.sqrt(np.mean(delta_param**2))
            if norm_delta_param < STOP_CONDITION_DELTA_PARAM:
                break


        # update parameters
        self.intrinsic['fu']      = param[0,0]
        self.intrinsic['fv']      = param[1,0]
        self.intrinsic['cu']      = param[2,0]
        self.intrinsic['cv']      = param[3,0]
        self.intrinsic['dist'][0] = param[4,0]
        self.intrinsic['dist'][1] = param[5,0]
        self.intrinsic['dist'][2] = param[6,0]

        for idx_image in range(num_images):
            extrin = param[num_intrinsics + idx_image*num_extrinsics:num_intrinsics + (idx_image+1)*num_extrinsics]
            qr = extrin[0,0]
            qx = extrin[1,0]
            qy = extrin[2,0]            
            qz = extrin[3,0]
            q_rot = Quaternion(np.array([qr, qx, qy, qz]))

            roll, pitch, yaw = q_rot.to_euler_angles()
            R = q_rot.to_rotation_matrix()
            rt = extrin[4:]
            t = np.matmul(R.T, rt)

            self.extrinsics[idx_image]['roll']  = roll
            self.extrinsics[idx_image]['pitch'] = pitch
            self.extrinsics[idx_image]['yaw']   = yaw
            self.extrinsics[idx_image]['tx']    = t[0,0]
            self.extrinsics[idx_image]['ty']    = t[1,0]
            self.extrinsics[idx_image]['tz']    = t[2,0]

    def compute_f(self, intrinsic, extrinsic, uv, xyz):
        """This method computes error vector.
        """
        u, v = uv
        x, y, z = xyz

        f_u = intrinsic[0]
        f_v = intrinsic[1]
        cu  = intrinsic[2]
        cv  = intrinsic[3]     
        k1  = intrinsic[4]
        k2  = intrinsic[5]
        k3  = intrinsic[6]
        
        qr  = extrinsic[0]
        qx  = extrinsic[1]
        qy  = extrinsic[2]
        qz  = extrinsic[3]
        t_x = extrinsic[4]
        t_y = extrinsic[5]
        t_z = extrinsic[6]

        t1 = (-2.0*qy**2 - 2.0*qz**2 + 1.0)
        t2 = (-2.0*qr*qz + 2.0*qx*qy)
        t3 = (2.0*qr*qy + 2.0*qx*qz)
        t4 = (-2.0*qr*qy + 2.0*qx*qz)
        t5 = (2.0*qr*qx + 2.0*qy*qz)
        t6 = (-2.0*qx**2 - 2.0*qy**2 + 1.0)
        t7 = (2.0*qr*qz + 2.0*qx*qy)
        t8 = (-2.0*qx**2 - 2.0*qz**2 + 1.0)
        t9 = (-2.0*qr*qx + 2.0*qy*qz)

        t101 = (t_x + x*t1 + y*t2 + z*t3)
        t102 = (t_z + x*t4 + y*t5 + z*t6)
        t103 = (t_y + x*t7 + y*t8 + z*t9)

        t1001 = (t101**2/t102**2 + t103**2/t102**2)

        e_u = cu + f_u*t101*(k1*t1001 + k2*t1001**2 + k3*t1001**3 + 1)/t102 - u
        e_v = cv + f_v*t103*(k1*t1001 + k2*t1001**2 + k3*t1001**3 + 1)/t102 - v
        e_q = qr**2 + qx**2 + qy**2 + qz**2 - 1

        f = np.array([[e_u], [e_v], [e_q]])

        return f

    def compute_J(self, intrinsic:np.ndarray, extrinsic:np.ndarray, xyz:np.ndarray):
        """This method computes Jacobian matrix
        """
        x, y, z = xyz

        f_u = intrinsic[0]
        f_v = intrinsic[1]

        k1  = intrinsic[4]
        k2  = intrinsic[5]
        k3  = intrinsic[6]
        
        qr  = extrinsic[0]
        qx  = extrinsic[1]
        qy  = extrinsic[2]
        qz  = extrinsic[3]
        t_x = extrinsic[4]
        t_y = extrinsic[5]
        t_z = extrinsic[6]
        
        t1 = (-2.0*qy**2 - 2.0*qz**2 + 1.0)
        t2 = (-2.0*qr*qz + 2.0*qx*qy)
        t3 = (2.0*qr*qy + 2.0*qx*qz)
        t4 = (-2.0*qr*qy + 2.0*qx*qz)
        t5 = (2.0*qr*qx + 2.0*qy*qz)
        t6 = (-2.0*qx**2 - 2.0*qy**2 + 1.0)
        t7 = (2.0*qr*qz + 2.0*qx*qy)
        t8 = (-2.0*qx**2 - 2.0*qz**2 + 1.0)
        t9 = (-2.0*qr*qx + 2.0*qy*qz)
        t10 = (-4.0*qx*y + 4.0*qy*x)
        t11 = (-4.0*qx*z + 4.0*qz*x)
        t12 = (4.0*qy*z - 4.0*qz*y)
        t13 = (4.0*qy*y + 4.0*qz*z)
        t14 = (-4.0*qr*y + 8.0*qx*z - 4.0*qz*x)
        t15 = (-4.0*qr*z - 8.0*qx*y + 4.0*qy*x)
        t16 = (4.0*qx*x + 4.0*qz*z)
        t17 = (4.0*qr*x + 8.0*qy*z - 4.0*qz*y)
        t18 = (4.0*qr*z + 4.0*qx*y - 8.0*qy*x)
        t19 = (-4.0*qx*x - 4.0*qy*y)
        t20 = (4.0*qr*x + 4.0*qy*z - 8.0*qz*y)
        t21 = (-4.0*qr*y + 4.0*qx*z - 8.0*qz*x)
        t22 = (2*t_x + 2*x*t1 + 2*y*t2 + 2*z*t3)
        t23 = (2*t_y + 2*x*t7 + 2*y*t8 + 2*z*t9)

        t101 = (t_x + x*t1 + y*t2 + z*t3)
        t102 = (t_z + x*t4 + y*t5 + z*t6)
        t103 = (t_y + x*t7 + y*t8 + z*t9)

        t1001 = (t101**2/t102**2 + t103**2/t102**2)
        
        t10001 = (k1*t1001 + k2*t1001**2 + k3*t1001**3 + 1)

        de_dfu = np.array([[t101*t10001/t102],
                           [0],
                           [0]])
        de_dfv = np.array([[0],
                           [t103*t10001/t102], 
                           [0]])
        de_dcu = np.array([[1], [0], [0]])
        de_dcv = np.array([[0], [1], [0]])
        de_dk1 = np.array([[f_u*t1001*t101/t102],
                           [f_v*t1001*t103/t102],
                           [0]])
        de_dk2 = np.array([[f_u*t1001**2*t101/t102],
                           [f_v*t1001**2*t103/t102],
                           [0]])
        de_dk3 = np.array([[f_u*t1001**3*t101/t102], 
                           [f_v*t1001**3*t103/t102],
                           [0]])

        de_dqr = np.array([[f_u*(-2.0*qx*y + 2.0*qy*x)*t101*t10001/t102**2 + f_u*(2.0*qy*z - 2.0*qz*y)*t10001/t102 + f_u*(k1*(t10*t101**2/t102**3 + t10*t103**2/t102**3 + t11*t103/t102**2 + t12*t101/t102**2) + k2*t1001*(2*t10*t101**2/t102**3 + 2*t10*t103**2/t102**3 + 2*t11*t103/t102**2 + 2*t12*t101/t102**2) + k3*t1001**2*(3*t10*t101**2/t102**3 + 3*t10*t103**2/t102**3 + 3*t11*t103/t102**2 + 3*t12*t101/t102**2))*t101/t102],
                           [f_v*(-2.0*qx*y + 2.0*qy*x)*t103*t10001/t102**2 + f_v*(-2.0*qx*z + 2.0*qz*x)*t10001/t102 + f_v*(k1*(t10*t101**2/t102**3 + t10*t103**2/t102**3 + t11*t103/t102**2 + t12*t101/t102**2) + k2*t1001*(2*t10*t101**2/t102**3 + 2*t10*t103**2/t102**3 + 2*t11*t103/t102**2 + 2*t12*t101/t102**2) + k3*t1001**2*(3*t10*t101**2/t102**3 + 3*t10*t103**2/t102**3 + 3*t11*t103/t102**2 + 3*t12*t101/t102**2))*t103/t102],
                           [2*qr]])
        de_dqx = np.array([[f_u*(2.0*qy*y + 2.0*qz*z)*t10001/t102 + f_u*(k1*(t13*t101/t102**2 + t14*t101**2/t102**3 + t14*t103**2/t102**3 + t15*t103/t102**2) + k2*t1001*(2*t13*t101/t102**2 + 2*t14*t101**2/t102**3 + 2*t14*t103**2/t102**3 + 2*t15*t103/t102**2) + k3*t1001**2*(3*t13*t101/t102**2 + 3*t14*t101**2/t102**3 + 3*t14*t103**2/t102**3 + 3*t15*t103/t102**2))*t101/t102 + f_u*(-2.0*qr*y + 4.0*qx*z - 2.0*qz*x)*t101*t10001/t102**2], 
                           [f_v*(k1*(t13*t101/t102**2 + t14*t101**2/t102**3 + t14*t103**2/t102**3 + t15*t103/t102**2) + k2*t1001*(2*t13*t101/t102**2 + 2*t14*t101**2/t102**3 + 2*t14*t103**2/t102**3 + 2*t15*t103/t102**2) + k3*t1001**2*(3*t13*t101/t102**2 + 3*t14*t101**2/t102**3 + 3*t14*t103**2/t102**3 + 3*t15*t103/t102**2))*t103/t102 + f_v*(-2.0*qr*y + 4.0*qx*z - 2.0*qz*x)*t103*t10001/t102**2 + f_v*(-2.0*qr*z - 4.0*qx*y + 2.0*qy*x)*t10001/t102],
                           [2*qx]])
        de_dqy = np.array([[f_u*(k1*(t16*t103/t102**2 + t17*t101**2/t102**3 + t17*t103**2/t102**3 + t18*t101/t102**2) + k2*t1001*(2*t16*t103/t102**2 + 2*t17*t101**2/t102**3 + 2*t17*t103**2/t102**3 + 2*t18*t101/t102**2) + k3*t1001**2*(3*t16*t103/t102**2 + 3*t17*t101**2/t102**3 + 3*t17*t103**2/t102**3 + 3*t18*t101/t102**2))*t101/t102 + f_u*(2.0*qr*x + 4.0*qy*z - 2.0*qz*y)*t101*t10001/t102**2 + f_u*(2.0*qr*z + 2.0*qx*y - 4.0*qy*x)*t10001/t102],
                           [f_v*(2.0*qx*x + 2.0*qz*z)*t10001/t102 + f_v*(k1*(t16*t103/t102**2 + t17*t101**2/t102**3 + t17*t103**2/t102**3 + t18*t101/t102**2) + k2*t1001*(2*t16*t103/t102**2 + 2*t17*t101**2/t102**3 + 2*t17*t103**2/t102**3 + 2*t18*t101/t102**2) + k3*t1001**2*(3*t16*t103/t102**2 + 3*t17*t101**2/t102**3 + 3*t17*t103**2/t102**3 + 3*t18*t101/t102**2))*t103/t102 + f_v*(2.0*qr*x + 4.0*qy*z - 2.0*qz*y)*t103*t10001/t102**2],
                           [2*qy]])
        de_dqz = np.array([[f_u*(-2.0*qx*x - 2.0*qy*y)*t101*t10001/t102**2 + f_u*(k1*(t19*t101**2/t102**3 + t19*t103**2/t102**3 + t20*t103/t102**2 + t21*t101/t102**2) + k2*t1001*(2*t19*t101**2/t102**3 + 2*t19*t103**2/t102**3 + 2*t20*t103/t102**2 + 2*t21*t101/t102**2) + k3*t1001**2*(3*t19*t101**2/t102**3 + 3*t19*t103**2/t102**3 + 3*t20*t103/t102**2 + 3*t21*t101/t102**2))*t101/t102 + f_u*(-2.0*qr*y + 2.0*qx*z - 4.0*qz*x)*t10001/t102],
                           [f_v*(-2.0*qx*x - 2.0*qy*y)*t103*t10001/t102**2 + f_v*(k1*(t19*t101**2/t102**3 + t19*t103**2/t102**3 + t20*t103/t102**2 + t21*t101/t102**2) + k2*t1001*(2*t19*t101**2/t102**3 + 2*t19*t103**2/t102**3 + 2*t20*t103/t102**2 + 2*t21*t101/t102**2) + k3*t1001**2*(3*t19*t101**2/t102**3 + 3*t19*t103**2/t102**3 + 3*t20*t103/t102**2 + 3*t21*t101/t102**2))*t103/t102 + f_v*(2.0*qr*x + 2.0*qy*z - 4.0*qz*y)*t10001/t102],
                           [2*qz]])
        de_dtx = np.array([[f_u*(k1*t22/t102**2 + 2*k2*t1001*t22/t102**2 + 3*k3*t1001**2*t22/t102**2)*t101/t102 + f_u*t10001/t102],
                           [f_v*(k1*t22/t102**2 + 2*k2*t1001*t22/t102**2 + 3*k3*t1001**2*t22/t102**2)*t103/t102],
                           [0]])
        de_dty = np.array([[f_u*(k1*t23/t102**2 + 2*k2*t1001*t23/t102**2 + 3*k3*t1001**2*t23/t102**2)*t101/t102],
                           [f_v*(k1*t23/t102**2 + 2*k2*t1001*t23/t102**2 + 3*k3*t1001**2*t23/t102**2)*t103/t102 + f_v*t10001/t102],
                           [0]])
        de_dtz = np.array([[f_u*(k1*(-2*t101**2/t102**3 - 2*t103**2/t102**3) + k2*(-4*t101**2/t102**3 - 4*t103**2/t102**3)*t1001 + k3*(-6*t101**2/t102**3 - 6*t103**2/t102**3)*t1001**2)*t101/t102 - f_u*t101*t10001/t102**2],
                           [f_v*(k1*(-2*t101**2/t102**3 - 2*t103**2/t102**3) + k2*(-4*t101**2/t102**3 - 4*t103**2/t102**3)*t1001 + k3*(-6*t101**2/t102**3 - 6*t103**2/t102**3)*t1001**2)*t103/t102 - f_v*t103*t10001/t102**2],
                           [0]])

        J_intrin = np.hstack((de_dfu, de_dfv, de_dcu, de_dcv, de_dk1, de_dk2, de_dk3))
        J_extrin = np.hstack((de_dqr, de_dqx, de_dqy, de_dqz, de_dtx, de_dty, de_dtz))

        return J_intrin, J_extrin


        