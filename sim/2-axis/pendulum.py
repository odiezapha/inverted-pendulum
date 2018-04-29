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

class Controller():
    def __init__(self, motor, wheel, pend):
        ## Config ##
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

    def get_system_state_matrix(self):
        """Return 9x9 system state matrix"""

        # Create empty matrix to describe system state
        A = zeros((9,9))

        # Q-dot terms
        A[0:3,3:6] = 0.5*identity(3)

        # X and Y gravity terms and Z rotation
        A[3:6,0:3] = matrix([
            [2*self.g*self.lp*self.Mp/(self.Ix), 0, 0],
            [0, 2*self.g*self.lp*self.Mp/(self.Iy), 0],
            [0, 0, 2*self.g*self.Fspin_z/(self.Iz)]
        ])
        # Friction in pendulum rotation on pendulum
        A[3:6,3:6] = matrix([
            [self.Bx/(self.Ix), 0, 0],
            [0, self.By/(self.Iy), 0],
            [0, 0, self.Bz/(self.Iz)]
        ])
        # Friction in wheel rotation on pendulum
        A[3:6,6:9] = matrix([
            [-self.Bw/(self.Ix), 0, 0],
            [0, -self.Bw/(self.Iy), 0],
            [0, 0, -self.Bw/(self.Iz)]
        ])
        # Friction in wheel rotation on wheel
        A[6:9,6:9] = -self.Bw/self.Iw*identity(3)

        return A

    def get_input_state_matrix(self):
        """Return 9x3 input state matrix"""

        # state matrix for control
        B = zeros((9,3))

        # control term impact on pendulum
        B[3:6,0:3] = matrix([
            [self.Km/(self.Ix)/self.Rm, 0, 0],
            [0, self.Km/(self.Iy)/self.Rm, 0],
            [0, 0, self.Km/(self.Iz)/self.Rm]
        ])
        # control term impact on wheel
        B[6:9,0:3] = self.Km/self.Iw/self.Rm*identity(3)

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

    def print_eigs(self, A, B, K):
        """Return Eigenvalues and Eigenvalues"""

        # compute the Eigenvalues and Eigenvectors
        eigVals, eigVecs = scipy.linalg.eig(A-B*K)
        print("Eigen Values: ")
        print(eigVals)
        print("Eigen Vectors: ")
        print(eigVecs)


class Sim():
    def __init__(self, motor, wheel, pend, K):
        self.k = K
        self.v = np.zeros(2) # initialize control voltages

        ## Config ##
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

    def deriv(self, u, t):
        """First order ODEs for integration"""

        # Generate control voltages, clip at saturation voltage
        Vx = np.clip(float(-self.k[0]*matrix((u)).T), 0, self.Vsat)
        Vy = np.clip(float(-self.k[1]*matrix((u)).T), 0, self.Vsat)
        Vz = 0

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
        x4_dt = self.Km*Vx/self.Rm/self.Ix + 2*self.g*self.lp*self.Mp/(self.Ix)*x1 + self.Bx/(self.Ix)*x4 - self.Bw/(self.Ix)*x7
        x5_dt = self.Km*Vy/self.Rm/self.Iy + 2*self.g*self.lp*self.Mp/(self.Iy)*x2 + self.By/(self.Iy)*x5 - self.Bw/(self.Iy)*x8
        x6_dt = self.Km*Vz/self.Rm/self.Iz + 2*self.g*self.Fspin_z/(self.Iz)*x3 + self.Bz/(self.Iz)*x6 - self.Bw/(self.Iz)*x9
        x7_dt = self.Km*Vx/self.Rm/self.Iw + -self.Bw/self.Iw*x7
        x8_dt = self.Km*Vy/self.Rm/self.Iw + -self.Bw/self.Iw*x8
        x9_dt = self.Km*Vz/self.Rm/self.Iw + -self.Bw/self.Iw*x9
        dxdt = [x1_dt, x2_dt, x3_dt, x4_dt, x5_dt, x6_dt, x7_dt, x8_dt, x9_dt]

        return dxdt

    def gen_data(self, system_init, dt, t_end):
        """ Use scipy.integrate.odeint to generate time series data"""

        x = system_init[:] # initialize system state
        t = np.arange(0.0, t_end, dt)

        ts_data = odeint(self.deriv, x, t)

        return ts_data, t

    def plot_timeseries(self, ts_data, t):
        plt.plot(t, ts_data[:, 0], 'b', label='q0(t)')
        plt.plot(t, ts_data[:, 1], 'g', label='q1(t)')
        plt.legend(loc='best')
        plt.xlabel('t')
        plt.grid()
        plt.show()

def main():

    # Set configuration
    motor = config_motor.maxon_70w
    wheel = config_wheel.estimate
    pend = config_pendulum.estimate

    C = Controller(motor, wheel, pend)

    A = C.get_system_state_matrix()
    B = C.get_input_state_matrix()
    R = C.get_cost_matrix(gains = [1, 1, 1])
    Q = C.get_effort_matrix(gains = [800, 800, 800, 800, 800, 800, 1, 1, 1])

    # Get LQR controller gains
    K = C.get_k_gain(A, B, Q, R)

    # C.print_eigs(A, B, K)

    # Set up simulation
    duration = 30 # sec
    time_step = 0.01 # sec
    sys_init = [
            0.03, # q0
            0.04, # q1
            0.,
            0.1, # w_x
            0.5, # w_y
            0.,
            0.,
            0.,
            0.]

    S = Sim(motor, wheel, pend, K)

    ts_data, t = S.gen_data(sys_init, time_step, duration)

    S.plot_timeseries(ts_data, t)

    return

if __name__== "__main__":
  main()
