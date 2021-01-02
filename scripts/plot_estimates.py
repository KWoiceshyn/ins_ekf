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
    x_cov = cov_array[:,4]
    y = states_array[:,5]
    y_cov = cov_array[:,5]
    z = states_array[:,6]
    z_cov = cov_array[:,6]

    Vx = states_array[:,7]
    Vx_cov = cov_array[:,7]
    Vy = states_array[:,8]
    Vy_cov = cov_array[:,8]
    Vz = states_array[:,9]
    Vz_cov = cov_array[:,9]

    Vb_normsq = Vx * Vx + Vy * Vy + Vz * Vz
    Vb = np.sqrt(Vb_normsq)

    plt.figure(1)
    plt.plot(x, y, c='r', label='path')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('2D Path')
    plt.axis('equal')
    plt.legend()

    plt.figure(2)
    plt.plot(t, Vb, c='r', label='velocity')
    plt.xlabel('sec')
    plt.ylabel('m/s')
    plt.title('Body Velocity')
    plt.legend()

    fig3, ax = plt.subplots(3)
    fig3.suptitle("XYZ Position")
    ax[0].plot(t, x, c='r', label='x_position')
    ax[0].plot(t, x + x_cov, c='r', ls='--', lw=0.5)
    ax[0].plot(t, x - x_cov, c='r', ls='--', lw=0.5)
    ax[0].set_ylabel('m')
    ax[0].legend()
    ax[1].plot(t, y, c='r', label='y_position')
    ax[1].plot(t, y + y_cov, c='r', ls='--', lw=0.5)
    ax[1].plot(t, y - y_cov, c='r', ls='--', lw=0.5)
    ax[1].set_ylabel('m')
    ax[1].legend()
    ax[2].plot(t, z, c='r', label='z_position')
    ax[2].plot(t, z + z_cov, c='r', ls='--', lw=0.5)
    ax[2].plot(t, z - z_cov, c='r', ls='--', lw=0.5)
    ax[2].set_ylabel('m')
    ax[2].legend()
    ax[2].set_xlabel('sec')

    fig4, ax2 = plt.subplots(3)
    fig4.suptitle("RPY")
    ax2[0].plot(t, roll, c='r', label='roll')
    ax2[0].plot(t, roll + roll_cov, c='r', ls='--', lw=0.5)
    ax2[0].plot(t, roll - roll_cov, c='r', ls='--', lw=0.5)
    ax2[0].set_ylabel('deg')
    ax2[0].legend()
    ax2[1].plot(t, pitch, c='r', label='pitch')
    ax2[1].plot(t, pitch + pitch_cov, c='r', ls='--', lw=0.5)
    ax2[1].plot(t, pitch - pitch_cov, c='r', ls='--', lw=0.5)
    ax2[1].set_ylabel('deg')
    ax2[1].legend()
    ax2[2].plot(t, yaw, c='r', label='yaw')
    ax2[2].plot(t, yaw + yaw_cov, c='r', ls='--', lw=0.5)
    ax2[2].plot(t, yaw - yaw_cov, c='r', ls='--', lw=0.5)
    ax2[2].set_ylabel('deg')
    ax2[2].legend()
    ax2[2].set_xlabel('sec')

    fig5, ax3 = plt.subplots(3)
    fig5.suptitle("Vxyz")
    ax3[0].plot(t, Vx, c='r', label='Vx')
    ax3[0].plot(t, Vx + Vx_cov, c='r', ls='--', lw=0.5)
    ax3[0].plot(t, Vx - Vx_cov, c='r', ls='--', lw=0.5)
    ax3[0].set_ylabel('m/s')
    ax3[0].legend()
    ax3[1].plot(t, Vy, c='r', label='Vy')
    ax3[1].plot(t, Vy + Vy_cov, c='r', ls='--', lw=0.5)
    ax3[1].plot(t, Vy - Vy_cov, c='r', ls='--', lw=0.5)
    ax3[1].set_ylabel('m/s')
    ax3[1].legend()
    ax3[2].plot(t, Vz, c='r', label='Vz')
    ax3[2].plot(t, Vz + Vz_cov, c='r', ls='--', lw=0.5)
    ax3[2].plot(t, Vz - Vz_cov, c='r', ls='--', lw=0.5)
    ax3[2].set_ylabel('m/s')
    ax3[2].legend()
    ax3[2].set_xlabel('sec')

    plt.show()


if __name__ == "__main__":
    main()
