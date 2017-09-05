from __future__ import division, print_function
from math import sin, cos, pi
from numpy import matrix, array
import scipy.linalg
import controlpy

R = [2.2]

Q = matrix([
    [10, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
])

m = 1        # mass of wheel+pendulum (kg)
Km = 0.07    # motor torque constant (Nm/A)
Rm = 6.69    # motor internal resistance (Ohms)
Vsat = 22    # motor max voltage
Iw = 0.03    # wheel MOI (kg*m^2)
Bw = 0.001   # friction in wheel bearing (Nms)
l = 0.33     # length of pendulum (m)
Ip = 0.05    # pendulum MOI (kg*m^2)
Bp = 0       # friction in pendulum bearing (Nms)
g = -9.81    # gravity (m/s^2)

A10 = g*l*m/(Iw+Ip) # should be multiplied by cos(x1d)
A11 = Bp/(Iw+Ip)
A13 = -Bw/(Iw+Ip)
A33 = -Bw/Iw

A = matrix([
    [0, 1, 0, 0],
    [A10, A11, 0, A13],
    [0, 0, 0, 1],
    [0, 0, 0, A33]
])
print("A: ")
print(A)

B1 = Km/((Iw+Ip)*Rm)
B3 = Km/(Iw*Rm)

B = matrix([
    [0],
    [B1],
    [0],
    [B3]
])
print("B: ")
print(B)

hurwitz_check = controlpy.analysis.is_hurwitz(A)

if hurwitz_check:
    print('Matrix A meets Hurwitz criteria')
else:
    print('Matrix A does not meet Hurwitz criteria')

uncontrollableModes = controlpy.analysis.uncontrollable_modes(A,B)

if not uncontrollableModes:
    print('System is controllable.')
else:
    print('System is uncontrollable. Uncontrollable modes are:')
    print(uncontrollableModes)

#solve the Ricatti equation
X = matrix(scipy.linalg.solve_continuous_are(A, B, Q, R))
print("Solution: ")
print(X)

#compute the LQR gain
K = matrix((R)*(B.T*X))
print("Gains: ")
print(K)

#compute the Eigenvalues and Eigenvectors
eigVals, eigVecs = scipy.linalg.eig(A-B*K)
print("Eigen Vectors: ")
print(eigVecs)
print("Eigen Values: ")
print(eigVals)


def average(x):
    x_i, k1, k2, k3, k4 = x
    return x_i + (k1 + 2.0*(k3 + k4) +  k2) / 6.0

def constrain(theta):
    theta = theta % (2*pi)
    if theta > pi:
        theta = -2*pi+theta
    return theta

def sat(Vsat, V):
    if abs(V) > Vsat:
        return Vsat * cmp(V, 0)
    return V

theta = []
class Pendulum(object):
    def __init__(self, dt, init_conds, end):
        self.dt = dt
        self.t = 0.0
        self.x = init_conds[:]
        self.v = 0
        self.end = end

    def derivative(self, u):
        V = sat(Vsat, self.control(u))
        #x1 = Theta Pendulum, x2 = W Pendulum, x3 = Theta Wheel, x4 = W Wheel
        x1, x2, x3, x4 = u
        x1_dt, x3_dt =  x2, x4
        x2_dt = (Km*V/Rm - Bw*x4 + Bp*x2 + g*l*m*sin(x1))/(Iw + Ip)
        x4_dt = (Km*V/Rm - Bw*x4)/(Iw)
        x = [x1_dt, x2_dt, x3_dt, x4_dt]
        return x

    def control(self, u):
        c = constrain(u[0])
        set_point = [0,0,0,0]
        if c>-pi/5 and c<pi/5:
            return float(-K*(matrix([c]+u[0:3])-matrix(set_point)).T)
        else:
            return 0

    def rk4_step(self, dt):
        dx = self.derivative(self.x)
        k2 = [ dx_i*dt for dx_i in dx ]

        xv = [x_i + delx0_i/2.0 for x_i, delx0_i in zip(self.x, k2)]
        k3 = [ dx_i*dt for dx_i in self.derivative(xv)]

        xv = [x_i + delx1_i/2.0 for x_i,delx1_i in zip(self.x, k3)]
        k4 = [ dx_i*dt for dx_i in self.derivative(xv) ]

        xv = [x_i + delx1_2 for x_i,delx1_2 in zip(self.x, k4)]
        k1 = [self.dt*i for i in self.derivative(xv)]

        self.v = sat(Vsat, self.control(self.x))
        self.t += dt
        self.x = map(average, zip(self.x, k1, k2, k3, k4))

    def integrate(self):
        x = []
        while self.t <= self.end:
            self.rk4_step(self.dt)
            x.append([self.t] + self.x + [self.v])
        return array(x)
