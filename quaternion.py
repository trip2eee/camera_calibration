"""
@fn     quaternion.py
@brief  Quaternion library (github.com/trip2eee/quatlib)
@author Jongmin Park
@date   July 18, 2022
"""
import numpy as np

class Quaternion:
    def __init__(self, q=None):

        if q is not None and isinstance(q, np.ndarray):
            self.q = q.copy()
        else:
            self.q = np.array([0] * 4, dtype=np.float32)
    
    def __getitem__(self, i:int):
        return self.q[i]

    def __str__(self):
        s = ""
        if abs(self.q[0]) > 0:
            s += "{}".format(self.q[0])
        
        if self.q[1] > 0:
            s += "+{}i".format(self.q[1])
        elif self.q[1] < 0:
            s += "{}i".format(self.q[1])

        if self.q[2] > 0:
            s += "+{}j".format(self.q[2])
        elif self.q[2] < 0:
            s += "{}j".format(self.q[2])

        if self.q[3] > 0:
            s += "+{}k".format(self.q[3])
        elif self.q[3] < 0:
            s += "{}k".format(self.q[3])

        return s

    def __add__(self, other:'Quaternion') -> 'Quaternion':
        q = self.q + other.q
        result = Quaternion(q)
        return result

    def __sub__(self, other:'Quaternion') -> 'Quaternion':
        q = self.q - other.q
        result = Quaternion(q)
        return result

    def __mul__(self, other:'Quaternion') -> 'Quaternion':
        a1 = self.q[0]
        b1 = self.q[1]
        c1 = self.q[2]
        d1 = self.q[3]

        a2 = other.q[0]
        b2 = other.q[1]
        c2 = other.q[2]
        d2 = other.q[3]

        a3 = (a1*a2) - (b1*b2) - (c1*c2) - (d1*d2) # real
        b3 = (a1*b2) + (b1*a2) + (c1*d2) - (d1*c2) # i
        c3 = (a1*c2) - (b1*d2) + (c1*a2) + (d1*b2) # j
        d3 = (a1*d2) + (b1*c2) - (c1*b2) + (d1*a2) # k

        q3 = np.array([a3, b3, c3, d3])
        result = Quaternion(q3)
        return result

    def conjugate(self) -> 'Quaternion':
        result = Quaternion()
        result.q[0] =  self.q[0]
        result.q[1] = -self.q[1]
        result.q[2] = -self.q[2]
        result.q[3] = -self.q[3]
        return result

    def norm(self):
        a = self.q[0]
        b = self.q[1]
        c = self.q[2]
        d = self.q[3]

        return np.sqrt(a**2 + b**2 + c**2 + d**2)

    def invert(self) -> 'Quaternion':
        result = Quaternion(self.q)
        norm = self.norm()

        result.q[0] /=  norm
        result.q[1] /= -norm
        result.q[2] /= -norm
        result.q[3] /= -norm
        
        return result
    
    def rotate(self, qr:'Quaternion') -> 'Quaternion':
        """This method rotates self by qr.
        """
        result = qr*self*qr.conjugate()
        return result
    
    def to_euler_angles(self):
        """This method computes Euler angles from the quaternion.
        """
        qr = self.q[0]
        qx = self.q[1]
        qy = self.q[2]
        qz = self.q[3]

        roll  = np.arctan2(2.0*((qr*qx) + (qy*qz)), 1.0 - 2.0*((qx*qx) + (qy*qy)))
        pitch =  np.arcsin(2.0*((qr*qy) - (qz*qx)))
        yaw   = np.arctan2(2.0*((qr*qz) + (qx*qy)), 1.0 - 2.0*((qy*qy) + (qz*qz)))

        return roll, pitch, yaw
    
    def to_quaternion(self, roll, pitch, yaw) -> 'Quaternion':
        """This method converts Euler angles to the quaternion.
        """
        c_roll  = np.cos(roll * 0.5)
        c_pitch = np.cos(pitch * 0.5)
        c_yaw   = np.cos(yaw * 0.5)

        s_roll  = np.sin(roll * 0.5)
        s_pitch = np.sin(pitch * 0.5)
        s_yaw   = np.sin(yaw * 0.5)

        qr = ((c_roll*c_pitch)*c_yaw) + ((s_roll*s_pitch)*s_yaw)
        qx = ((s_roll*c_pitch)*c_yaw) - ((c_roll*s_pitch)*s_yaw)
        qy = ((c_roll*s_pitch)*c_yaw) + ((s_roll*c_pitch)*s_yaw)
        qz = ((c_roll*c_pitch)*s_yaw) - ((s_roll*s_pitch)*c_yaw)
        
        self.q[0] = qr
        self.q[1] = qx
        self.q[2] = qy
        self.q[3] = qz

        return self

    def to_rotation_matrix(self) -> np.ndarray:
        """This method computes 3x3 rotation matrix from the quaternion.
        """
        qr = self.q[0]
        qx = self.q[1]
        qy = self.q[2]
        qz = self.q[3]

        r0 = 1.0 - 2.0*(qy*qy + qz*qz)
        r1 = 2.0*(qx*qy - qr*qz)
        r2 = 2.0*(qr*qy + qx*qz)

        r3 = 2.0*(qx*qy + qr*qz)
        r4 = 1.0 - 2.0*(qx*qx + qz*qz)
        r5 = 2.0*(qy*qz - qr*qx)

        r6 = 2.0*(qx*qz - qr*qy)
        r7 = 2.0*(qr*qx + qy*qz)
        r8 = 1.0 - 2.0*(qx*qx + qy*qy)

        return np.array([[r0, r1, r2],[r3, r4, r5], [r6, r7, r8]])


    def rotation_axis(self):
        """This method computes rotation axis of the quaternion.
        """
        qx = self.q[1]
        qy = self.q[2]
        qz = self.q[3]

        nv = np.sqrt(qx**2 + qy**2 + qz**2)

        ax = self.q[1] / nv
        ay = self.q[2] / nv
        az = self.q[3] / nv

        return ax, ay, az
    
    def rotation_angle(self):
        """This method computes rotation angle of the quaternion.
        """
        qr = self.q[0]
        qx = self.q[1]
        qy = self.q[2]
        qz = self.q[3]

        nv = np.sqrt(qx**2 + qy**2 + qz**2)
        return 2.0*np.arctan2(nv, qr)

