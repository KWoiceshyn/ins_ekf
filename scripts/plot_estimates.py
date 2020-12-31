#! usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def main():

    states_file = open("../states_out.csv", "rb")
    states_array = np.loadtxt(states_file, delimiter=",")

    r2d = 180 / np.pi
    t = states_array[:,0]
    x = states_array[:,1]
    y = states_array[:,2]
    roll = r2d * states_array[:,3]
    pitch = r2d * states_array[:,4]
    yaw = r2d * states_array[:,5]

    Vx = states_array[:,6]
    Vy = states_array[:,7]
    Vz = states_array[:,8]

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
    plt.xlabel('sec')
    plt.ylabel('deg')
    plt.title('Roll')
    plt.legend()
    
    plt.figure(3)
    
    plt.plot(t, pitch, c='r', label='pitch')
    plt.xlabel('sec')
    plt.ylabel('deg')
    plt.title('Pitch')
    plt.legend()

    plt.figure(4)

    plt.plot(t, yaw, c='r', label='yaw')
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
