from __future__ import division, print_function
import math
from math import sin, cos, pi
import numpy as np
from numpy import matrix, array, zeros, identity
import scipy.linalg
from scipy.integrate import odeint
from matplotlib import pyplot as plt

import config_motor
import config_wheel
import config_pendulum

class Pendulum:
    def __init__(self, motor, wheel, pend):
        self.Km = motor['Km']   # motor torque constant (Nm/A)
        self.Rm = motor['Rm']   # motor internal resistance (Ohms)
        self.Vsat = motor['Vsat']  # motor max voltage
        self.Iw = wheel['Iw']   # wheel MOI (kg*m^2)
        self.Bw = wheel['Bw']   # friction in wheel bearing (Nms)
        self.Mr = pend['Mr']    # mass of rod (kg)
        self.Yr = pend['Yr']    # normalized distance of rod along rod
        self.Mw = pend['Mw']    # mass of wheel (kg)
        self.Yw = pend['Yw']    # normalized distance of wheel along rod
        self.lp = pend['lp']    # length of pendulum (m)
        self.Ix = pend['Ix']    # pendulum MOI (kg*m^2)
        self.Bx = pend['Bx']    # friction in pendulum bearing (Nms)
        self.Iy = pend['Iy']    # pendulum MOI (kg*m^2)
        self.By = pend['By']    # friction in pendulum bearing (Nms)
        self.Iz = pend['Iz']    # pendulum MOI (kg*m^2)
        self.Bz = pend['Bz']    # friction in pendulum bearing (Nms)
        self.Fspin_z = pend['Fspin_z']  # force spinning around rod axis, term needed for stability

        self.Mp = self.Mw*self.Yw + self.Mr*self.Yr # pendulum's weighted mass

        self.g = -9.81    # gravity (m/s^2)

class Controller():
    def __init__(self, Pendulum):
        self.params = Pendulum.__dict__

    def get_system_state_matrix(self):
        """Return 9x9 system state matrix"""

        # Create empty matrix to describe system state
        A = zeros((9,9))

        # Q-dot terms
        A[0:3,3:6] = 0.5*identity(3)

        # X and Y gravity terms and Z rotation
        A[3:6,0:3] = matrix([
            [2*self.params['g']*self.params['lp']*self.params['Mp']/(self.params['Ix']), 0, 0],
            [0, 2*self.params['g']*self.params['lp']*self.params['Mp']/(self.params['Iy']), 0],
            [0, 0, 2*self.params['g']*self.params['Fspin_z']/(self.params['Iz'])]
        ])
        # Friction in pendulum rotation on pendulum
        A[3:6,3:6] = matrix([
            [self.params['Bx']/(self.params['Ix']), 0, 0],
            [0, self.params['By']/(self.params['Iy']), 0],
            [0, 0, self.params['Bz']/(self.params['Iz'])]
        ])
        # Friction in wheel rotation on pendulum
        A[3:6,6:9] = matrix([
            [-self.params['Bw']/(self.params['Ix']), 0, 0],
            [0, -self.params['Bw']/(self.params['Iy']), 0],
            [0, 0, -self.params['Bw']/(self.params['Iz'])]
        ])
        # Friction in wheel rotation on wheel
        A[6:9,6:9] = -self.params['Bw']/self.params['Iw']*identity(3)

        return A

    def get_input_state_matrix(self):
        """Return 9x3 input state matrix"""

        # state matrix for control
        B = zeros((9,3))

        # control term impact on pendulum
        B[3:6,0:3] = matrix([
            [self.params['Km']/(self.params['Ix'])/self.params['Rm'], 0, 0],
            [0, self.params['Km']/(self.params['Iy'])/self.params['Rm'], 0],
            [0, 0, self.params['Km']/(self.params['Iz'])/self.params['Rm']]
        ])
        # control term impact on wheel
        B[6:9,0:3] = self.params['Km']/self.params['Iw']/self.params['Rm']*identity(3)

        return B

    def get_cost_matrix(self, gains = np.ones(3)):
        """Return 3x3 cost optization matrix"""

        # gain matrix for cost optimization
        R = np.diag(gains)

        return R

    def get_effort_matrix(self, gains = 500*np.ones(9)):
        """Return 9x9 effort optimization matrix"""

        # gain matrix for effort optimization
        Q = np.diag(gains)

        return Q

    def get_k_gain(self, A, B, Q, R, output = 0):
        """Return K gain matrix"""

        # solve the Ricatti equation
        X = matrix(scipy.linalg.solve_continuous_are(A, B, Q, R))

        # compute the LQR gain
        K = matrix((R)*(B.T*X))

        if output:
            print("X-axis K Matrix Gains: ")
            print(K[0])
            print("Y-axis K Matrix Gains: ")
            print(K[1])

        return K

    def print_eigs(self, A, B, K, vals = 1, vecs = 0):
        """Return Eigenvalues and Eigenvalues of controllability matrix"""

        # compute the Eigenvalues and Eigenvectors
        eigVals, eigVecs = scipy.linalg.eig(A-B*K)

        if vals == 1:
            print("Eigen Values: ")
            print(eigVals)
        if vecs == 1:
            print("Eigen Vectors: ")
            print(eigVecs)


class Sim:
    def __init__(self, Pendulum):
        self.params = Pendulum.__dict__
        self.v = np.zeros(3) # initialize control voltages

    def deriv(self, u, t, K):
        """First order ODEs for integration"""

        # Generate control voltages, clip at saturation voltage
        self.v[0] = np.clip(float(-K[0]*matrix((u)).T), 0, self.params['Vsat'])
        self.v[1] = np.clip(float(-K[1]*matrix((u)).T), 0, self.params['Vsat'])
        self.v[2] = 0

        '''
        Definitions:
        x1 = q0
        x2 = q1
        x3 = q2
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
        x4_dt = self.params['Km']*self.v[0]/self.params['Rm']/self.params['Ix'] \
            + 2*self.params['g']*self.params['lp']*self.params['Mp']/(self.params['Ix'])*x1 \
            + self.params['Bx']/(self.params['Ix'])*x4 - self.params['Bw']/(self.params['Ix'])*x7
        x5_dt = self.params['Km']*self.v[1]/self.params['Rm']/self.params['Iy'] \
            + 2*self.params['g']*self.params['lp']*self.params['Mp']/(self.params['Iy'])*x2 \
            + self.params['By']/(self.params['Iy'])*x5 - self.params['Bw']/(self.params['Iy'])*x8
        x6_dt = self.params['Km']*self.v[2]/self.params['Rm']/self.params['Iz'] \
            + 2*self.params['g']*self.params['Fspin_z']/(self.params['Iz'])*x3 \
            + self.params['Bz']/(self.params['Iz'])*x6 - self.params['Bw']/(self.params['Iz'])*x9
        x7_dt = self.params['Km']*self.v[0]/self.params['Rm']/self.params['Iw'] \
            + -self.params['Bw']/self.params['Iw']*x7
        x8_dt = self.params['Km']*self.v[1]/self.params['Rm']/self.params['Iw'] \
            + -self.params['Bw']/self.params['Iw']*x8
        x9_dt = self.params['Km']*self.v[2]/self.params['Rm']/self.params['Iw'] \
            + -self.params['Bw']/self.params['Iw']*x9

        dxdt = [x1_dt, x2_dt, x3_dt, x4_dt, x5_dt, x6_dt, x7_dt, x8_dt, x9_dt]

        return dxdt

    def gen_data(self, system_init, dt, t_end, K):
        """ Use scipy.integrate.odeint to generate time series data"""

        x = system_init[:] # initialize system state
        t = np.arange(0.0, t_end, dt)

        ts_data = odeint(self.deriv, x, t, args=(K,)) # extra arguments need to be a tuple

        return ts_data, t

    def plot_timeseries(self, ts_data, t):
        plt.plot(t, ts_data[:, 0], 'b', label='q0(t)')
        plt.plot(t, ts_data[:, 1], 'g', label='q1(t)')
        plt.legend(loc='best')
        plt.xlabel('t (sec)')
        plt.ylabel('pendulum angle')
        plt.grid()
        plt.show()

    def plot_parametric(self, ts_data):
        plt.plot(ts_data[:, 3], ts_data[:, 0], 'b', label='q0, wheel-x')
        plt.plot(ts_data[:, 4], ts_data[:, 1], 'g', label='q1, wheel-y')
        plt.legend(loc='best')
        plt.xlabel('wheel speed (rad/sec)')
        plt.ylabel('pendulum angle')
        plt.grid()
        plt.show()

def main():

    # Set configuration
    motor = config_motor.maxon_70w
    wheel = config_wheel.estimate
    pend = config_pendulum.estimate

    # Create pendulum object
    pend = Pendulum(motor, wheel, pend)

    # Create controller based on pendulum
    C = Controller(pend)

    A = C.get_system_state_matrix()
    B = C.get_input_state_matrix()
    R = C.get_cost_matrix(gains = [1, 1, 1])
    Q = C.get_effort_matrix(gains = [800, 800, 800, 800, 800, 800, 1, 1, 1])

    # Get LQR controller gains
    K = C.get_k_gain(A, B, Q, R, output=1)

    # Set up simulation
    duration = 15 # sec
    time_step = 0.01 # sec
    sys_init = [
            -0.08, # q0
            0.04, # q1
            0.,
            0.1, # w_x
            -0.3, # w_y
            0.,
            0.,
            0.,
            0.]

    # Create sim of pendulum
    S = Sim(pend)

    ts_data, t = S.gen_data(sys_init, time_step, duration, K)

    # Plot data
    S.plot_timeseries(ts_data, t)
    S.plot_parametric(ts_data)

    return

if __name__== "__main__":
  main()
