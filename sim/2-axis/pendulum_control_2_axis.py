from __future__ import division, print_function
import math
from math import sin, cos, pi
import numpy as np
from numpy import matrix, array, zeros, identity
import scipy.linalg
import config_motor
import config_wheel
import config_pendulum

motor = config_motor.maxon_70w
wheel = config_wheel.estimate
pend = config_pendulum.estimate

## Config ##
Km = motor['Km']   # motor torque constant (Nm/A)
Rm = motor['Rm']   # motor internal resistance (Ohms)
Vsat = motor['Vsat']  # motor max voltage
Iw = wheel['Iw']   # wheel MOI (kg*m^2)
Bw = wheel['Bw']   # friction in wheel bearing (Nms)
Mr = pend['Mr']    # mass of rod (kg)
Yr = pend['Yr']    # normalized distance of rod along rod
Mw = pend['Mw']    # mass of wheel (kg)
Yw = pend['Yw']    # normalized distance of wheel along rod
lp = pend['lp']    # length of pendulum (m)
Ix = pend['Ix']    # pendulum MOI (kg*m^2)
Bx = pend['Bx']    # friction in pendulum bearing (Nms)
Iy = pend['Iy']    # pendulum MOI (kg*m^2)
By = pend['By']    # friction in pendulum bearing (Nms)
Iz = pend['Iz']    # pendulum MOI (kg*m^2)
Bz = pend['Bz']    # friction in pendulum bearing (Nms)
Fspin_z = pend['Fspin_z']  # force spinning around rod axis, term needed for stability

Mp = Mw*Yw + Mr*Yr #pendulum's weighted mass
g = -9.81    # gravity (m/s^2)

#state matrix for motion
A = zeros((9,9))

#quaternion derivative terms
A[0:3,3:6] = 0.5*identity(3)
#X and Y gravity terms and Z rotation
A[3:6,0:3] = matrix([
    [2*g*lp*Mp/(Ix), 0, 0],
    [0, 2*g*lp*Mp/(Iy), 0],
    [0, 0, 2*g*Fspin_z/(Iz)]
])
#friction in pendulum rotation on pendulum
A[3:6,3:6] = matrix([
    [Bx/(Ix), 0, 0],
    [0, By/(Iy), 0],
    [0, 0, Bz/(Iz)]
])
#friction in wheel rotation on pendulum
A[3:6,6:9] = matrix([
    [-Bw/(Ix), 0, 0],
    [0, -Bw/(Iy), 0],
    [0, 0, -Bw/(Iz)]
])
#friction in wheel rotation on wheel
A[6:9,6:9] = -Bw/Iw*identity(3)

#state matrix for control
B = zeros((9,3))

#control term impact on pendulum
B[3:6,0:3] = matrix([
    [Km/(Ix)/Rm, 0, 0],
    [0, Km/(Iy)/Rm, 0],
    [0, 0, Km/(Iz)/Rm]
])
#control term impact on wheel
B[6:9,0:3] = Km/Iw/Rm*identity(3)

# gain matrix for cost optimization
R = matrix([
    [1,0,0],
    [0,1,0],
    [0,0,1]
])

# gain matrix for effort optimization
Q = 500*identity(9)

#solve the Ricatti equation
X = matrix(scipy.linalg.solve_continuous_are(A, B, Q, R))

#compute the LQR gain
K = matrix((R)*(B.T*X))
print("X Gains: ")
print(K[0])
print("Y Gains: ")
print(K[1])

#compute the Eigenvalues and Eigenvectors
eigVals, eigVecs = scipy.linalg.eig(A-B*K)
print("Eigen Values: ")
print(eigVals)

def average(x):
    x_i, k1, k2, k3, k4 = x
    return x_i + (k1 + 2.0*(k3 + k4) +  k2) / 6.0

def sat(Vsat, V):
    if abs(V) > Vsat:
        return Vsat * cmp(V, 0)
    return V

class Pendulum(object):
    def __init__(self, dt, init_conds, end):
        self.dt = dt
        self.t = 0.0
        self.x = init_conds[:]
        self.v_x = 0
        self.v_y = 0
        self.X = 0
        self.Y = 0
        self.Z = 0
        self.end = end

    def derivative(self, u):
        Vx = sat(Vsat, self.control_x(u))
        Vy = sat(Vsat, self.control_y(u))
        Vz = 0
        '''
        x1 = q1
        x2 = q2
        x3 = q3
        x4 = W rod X
        x5 = W rod Y
        x6 = W rod Z (not in use)
        x7 = W wheel X
        x8 = W wheel Y
        x9 = W wheel Z (not in use)
        '''
        x1, x2, x3, x4, x5, x6, x7, x8, x9 = u
        x1_dt = 0.5*x4
        x2_dt = 0.5*x5
        x3_dt = 0.5*x6
        x4_dt = Km*Vx/Rm/Ix + 2*g*lp*Mp/(Ix)*x1 + Bx/(Ix)*x4 - Bw/(Ix)*x7
        x5_dt = Km*Vy/Rm/Iy + 2*g*lp*Mp/(Iy)*x2 + By/(Iy)*x5 - Bw/(Iy)*x8
        x6_dt = Km*Vz/Rm/Iz + 2*g*Fspin_z/(Iz)*x3 + Bz/(Iz)*x6 - Bw/(Iz)*x9
        x7_dt = Km*Vx/Rm/Iw + -Bw/Iw*x7
        x8_dt = Km*Vy/Rm/Iw + -Bw/Iw*x8
        x9_dt = Km*Vz/Rm/Iw + -Bw/Iw*x9
        x = [x1_dt, x2_dt, x3_dt, x4_dt, x5_dt, x6_dt, x7_dt, x8_dt, x9_dt]
        return x

    def control_x(self, u):
        #use the x-axis K values and inputs, u, to calculuate control
        return float(-K[0]*matrix((u)).T)

    def control_y(self, u):
        #use the y-axis K values and inputs, u, to calculuate control
        return float(-K[1]*(matrix(u)).T)

    def rk4_step(self, dt):
        dx = self.derivative(self.x)
        k2 = [ dx_i*dt for dx_i in dx ]

        xv = [x_i + delx0_i/2.0 for x_i, delx0_i in zip(self.x, k2)]
        k3 = [ dx_i*dt for dx_i in self.derivative(xv)]

        xv = [x_i + delx1_i/2.0 for x_i,delx1_i in zip(self.x, k3)]
        k4 = [ dx_i*dt for dx_i in self.derivative(xv) ]

        xv = [x_i + delx1_2 for x_i,delx1_2 in zip(self.x, k4)]
        k1 = [self.dt*i for i in self.derivative(xv)]

        self.v_x = sat(Vsat, self.control_x(self.x))
        self.v_y = sat(Vsat, self.control_y(self.x))
        self.t += dt
        self.x = map(average, zip(self.x, k1, k2, k3, k4))

    def q2euler(self, x, y, z, w):
	ysqr = y*y
	t0 = +2.0 * (w * x + y*z)
	t1 = +1.0 - 2.0 * (x*x + ysqr)
	self.X = math.degrees(math.atan2(t0, t1))
	t2 = +2.0 * (w*y - z*x)
	t2 =  1 if t2 > 1 else t2
	t2 = -1 if t2 < -1 else t2
	self.Y = math.degrees(math.asin(t2))
	t3 = +2.0 * (w * z + x*y)
	t4 = +1.0 - 2.0 * (ysqr + z*z)
	self.Z = math.degrees(math.atan2(t3, t4))

    def integrate(self):
        x = []
        while self.t <= self.end:
            self.rk4_step(self.dt)
            self.q2euler(self.x[0], self.x[1], self.x[2], math.sqrt(1-np.power(self.x[0],2)-np.power(self.x[1],2)-np.power(self.x[2],2)))
            x.append([self.t] + self.x + [self.v_x] + [self.v_y] + [self.X] + [self.Y] + [self.Z])
        return array(x)
