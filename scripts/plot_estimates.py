#! usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def main():

    states_file = open("../states_out.csv", "rb")
    cov_file = open("../cov_out.csv", "rb")
    states_array = np.loadtxt(states_file, delimiter=",")
    cov_array = np.loadtxt(cov_file, delimiter=",")

    r2d = 180 / np.pi
    t = states_array[:,0]
    roll = r2d * states_array[:,1]
    roll_cov = r2d * cov_array[:,1]
    pitch = r2d * states_array[:,2]
    pitch_cov = r2d * cov_array[:,2]
    yaw = r2d * states_array[:,3]
    yaw_cov = r2d * cov_array[:,3]

    x = states_array[:,4]
    y = states_array[:,5]

    Vx = states_array[:,7]
    Vy = states_array[:,8]
    Vz = states_array[:,9]

    Vb_normsq = Vx * Vx + Vy * Vy + Vz * Vz
    Vb = np.sqrt(Vb_normsq)

    plt.figure(1)

    plt.plot(x, y, c='r', label='path')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Position')
    plt.axis('equal')
    plt.legend()

    plt.figure(2)
    
    plt.plot(t, roll, c='r', label='roll')
    plt.plot(t, roll + roll_cov, c='r', ls='--', lw=0.5)
    plt.plot(t, roll - roll_cov, c='r', ls='--', lw=0.5)
    plt.xlabel('sec')
    plt.ylabel('deg')
    plt.title('Roll')
    plt.legend()
    
    plt.figure(3)
    
    plt.plot(t, pitch, c='r', label='pitch')
    plt.plot(t, pitch + pitch_cov, c='r', ls='--', lw=0.5)
    plt.plot(t, pitch - pitch_cov, c='r', ls='--', lw=0.5)
    plt.xlabel('sec')
    plt.ylabel('deg')
    plt.title('Pitch')
    plt.legend()

    plt.figure(4)

    plt.plot(t, yaw, c='r', label='yaw')
    plt.plot(t, yaw + yaw_cov, c='r', ls='--', lw=0.5)
    plt.plot(t, yaw - yaw_cov, c='r', ls='--', lw=0.5)
    plt.xlabel('sec')
    plt.ylabel('deg')
    plt.title('Yaw')
    plt.legend()

    plt.figure(5)

    plt.plot(t, Vb, c='r', label='velocity')
    plt.xlabel('sec')
    plt.ylabel('m/s')
    plt.title('Body Velocity')
    plt.legend()

    plt.show()


if __name__ == "__main__":
    main()
