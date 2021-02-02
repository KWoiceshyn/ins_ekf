#! usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt


def plot_data(states_filename, cov_filename, meas_filename):

    states_file = open(states_filename, "rb")
    cov_file = open(cov_filename, "rb")
    meas_file = open(meas_filename)
    states_array = np.loadtxt(states_file, delimiter=",")
    cov_array = np.loadtxt(cov_file, delimiter=",")
    meas_array = np.loadtxt(meas_file, delimiter=",")

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

    bGp = r2d * states_array[:,10]
    bGp_cov = r2d * cov_array[:,10]
    bGq = r2d * states_array[:,11]
    bGq_cov = r2d * cov_array[:,11]
    bGr = r2d * states_array[:,12]
    bGr_cov = r2d * cov_array[:,12]

    bAx = states_array[:,13]
    bAx_cov = cov_array[:,13]
    bAy = states_array[:,14]
    bAy_cov = cov_array[:,14]
    bAz = states_array[:,15]
    bAz_cov = cov_array[:,15]

    Vb_normsq = Vx * Vx + Vy * Vy + Vz * Vz
    Vb = np.sqrt(Vb_normsq) # body velocity is norm of inertial frame velocity

    p_meas = r2d * meas_array[:,0]
    q_meas = r2d * meas_array[:,1]
    r_meas = r2d * meas_array[:,2]
    ax_meas = meas_array[:,3]
    ay_meas = meas_array[:,4]
    az_meas = meas_array[:,5]

    x_meas = meas_array[:,6]
    y_meas = meas_array[:,7]
    z_meas = meas_array[:,8]
    v_meas = meas_array[:,9]
    yaw_meas = r2d * meas_array[:,10]

    plt.figure(1)
    plt.plot(x, y, c='r', label='path')
    plt.plot(x_meas, y_meas, c='g', label='gps')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('2D Path')
    plt.axis('equal')
    plt.legend()

    plt.figure(2)
    plt.plot(t, Vb, c='r', label='velocity')
    plt.plot(t, v_meas, c='g', label='gps')
    plt.xlabel('sec')
    plt.ylabel('m/s')
    plt.title('Body Velocity')
    plt.legend()

    fig3, ax = plt.subplots(3)
    fig3.suptitle("XYZ Position")
    ax[0].plot(t, x, c='r', label='x_position')
    ax[0].plot(t, x + x_cov, c='r', ls='--', lw=0.5)
    ax[0].plot(t, x - x_cov, c='r', ls='--', lw=0.5)
    ax[0].plot(t, x_meas, c='g', label='gps')
    ax[0].set_ylabel('m')
    ax[0].legend()
    ax[1].plot(t, y, c='r', label='y_position')
    ax[1].plot(t, y + y_cov, c='r', ls='--', lw=0.5)
    ax[1].plot(t, y - y_cov, c='r', ls='--', lw=0.5)
    ax[1].plot(t, y_meas, c='g', label='gps')
    ax[1].set_ylabel('m')
    ax[1].legend()
    ax[2].plot(t, z, c='r', label='z_position')
    ax[2].plot(t, z + z_cov, c='r', ls='--', lw=0.5)
    ax[2].plot(t, z - z_cov, c='r', ls='--', lw=0.5)
    ax[2].plot(t, z_meas, c='g', label='gps')
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
    ax2[2].plot(t, yaw_meas, c='g', label='gps')
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

    fig6, ax4 = plt.subplots(3)
    fig6.suptitle("biasPQR")
    ax4[0].plot(t, bGp, c='r', label='bP')
    ax4[0].plot(t, bGp + bGp_cov, c='r', ls='--', lw=0.5)
    ax4[0].plot(t, bGp - bGp_cov, c='r', ls='--', lw=0.5)
    ax4[0].set_ylabel('deg/s')
    ax4[0].legend()
    ax4[1].plot(t, bGq, c='r', label='bQ')
    ax4[1].plot(t, bGq + bGq_cov, c='r', ls='--', lw=0.5)
    ax4[1].plot(t, bGq - bGq_cov, c='r', ls='--', lw=0.5)
    ax4[1].set_ylabel('deg/s')
    ax4[1].legend()
    ax4[2].plot(t, bGr, c='r', label='bR')
    ax4[2].plot(t, bGr + bGr_cov, c='r', ls='--', lw=0.5)
    ax4[2].plot(t, bGr - bGr_cov, c='r', ls='--', lw=0.5)
    ax4[2].set_ylabel('deg/s')
    ax4[2].legend()
    ax4[2].set_xlabel('sec')

    fig7, ax5 = plt.subplots(3)
    fig7.suptitle("biasAxyz")
    ax5[0].plot(t, bAx, c='r', label='bAx')
    ax5[0].plot(t, bAx + bAx_cov, c='r', ls='--', lw=0.5)
    ax5[0].plot(t, bAx - bAx_cov, c='r', ls='--', lw=0.5)
    #ax5[0].plot(t, ax_meas, c='g', label='Ax')
    ax5[0].set_ylabel('m/s/s')
    ax5[0].legend()
    ax5[1].plot(t, bAy, c='r', label='bAy')
    ax5[1].plot(t, bAy + bAy_cov, c='r', ls='--', lw=0.5)
    ax5[1].plot(t, bAy - bAy_cov, c='r', ls='--', lw=0.5)
    #ax5[1].plot(t, ay_meas, c='g', label='Az')
    ax5[1].set_ylabel('m/s/s')
    ax5[1].legend()
    ax5[2].plot(t, bAz, c='r', label='bAz')
    ax5[2].plot(t, bAz + bAz_cov, c='r', ls='--', lw=0.5)
    ax5[2].plot(t, bAz - bAz_cov, c='r', ls='--', lw=0.5)
    #ax5[2].plot(t, az_meas, c='g', label='Az')
    ax5[2].set_ylabel('m/s/s')
    ax5[2].legend()
    ax5[2].set_xlabel('sec')

    fig8, ax6 = plt.subplots(3)
    fig8.suptitle("measPQR")
    ax6[0].plot(t, p_meas, c='g', label='p')
    ax6[0].set_ylabel('deg/s')
    ax6[0].legend()
    ax6[1].plot(t, q_meas, c='g', label='q')
    ax6[1].set_ylabel('deg/s')
    ax6[1].legend()
    ax6[2].plot(t, r_meas, c='g', label='r')
    ax6[2].set_ylabel('deg/s')
    ax6[2].legend()
    ax6[2].set_xlabel('sec')

    fig9, ax7 = plt.subplots(3)
    fig9.suptitle("measAxyz")
    ax7[0].plot(t, ax_meas, c='g', label='Ax')
    ax7[0].set_ylabel('m/s/s')
    ax7[0].legend()
    ax7[1].plot(t, ay_meas, c='g', label='Ay')
    ax7[1].set_ylabel('m/s/s')
    ax7[1].legend()
    ax7[2].plot(t, az_meas, c='g', label='Az')
    ax7[2].set_ylabel('m/s/s')
    ax7[2].legend()
    ax7[2].set_xlabel('sec')

    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Load files for plotting')
    parser.add_argument('-s', '--states_file', type=str, help='.csv file with EKF state estimates')
    parser.add_argument('-c', '--cov_file', type=str, help='.csv file with EKF state covariances')
    parser.add_argument('-m', '--meas_file', type=str, help='.csv file with measurement inputs to EKF')
    args = parser.parse_args()

    plot_data(args.states_file, args.cov_file, args.meas_file)

if __name__ == "__main__":
    main()
