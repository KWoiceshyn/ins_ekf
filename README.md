# INS EKF
I wrote this to experiment with the idea of an error-state EKF, to estimate the state of a vehicle that has an IMU and a GPS sensor. These sensors are assumed to be fixed to the vehicle body at an offset to the vehicle's center of rotation, and a relative orientation, which can be specified by (x,y,z) vectors and quaternions, respectively.

The IMU provides 6 measurement channels, consisting of 3-axis acceleration and 3-axis angular rates. The GPS provides 5 measurement channels, which consist of (x,y,z) position, velocity over ground, and track over ground (heading).

The EKF estimates the following states of the vehicle:
- Orientation (roll, pitch, yaw)
- Position (x,y,z)
- Inertial frame velocity (x,y,z)
- Gyro biases (x3)
- Accelerometer biases (x3)

The IMU measurements are used in the prediction step. The GPS measurements are used in the correction step, which occurs at a slower rate. 

## Usage
The input is a .txt or .csv file with measurement values, in the following columns:
| roll_rate | pitch_rate | yaw_rate | accel_x | accel_y | accel_z | x_pos | y_pos | z_pos | velocity | heading |
| -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- |


The first 6 measurements are from the IMU, and the last 5 are from the GPS.

For exampe, to run the EKF with a .csv file as input:
```
./main --mf measurements_file.csv
```

There is a script `plot_estimates.py` to plot the state estimates and covariances, used as follows:
```
python3 plot_estimates.py -s ../states_out.csv -c ../cov_out.csv -m ../measurements_file.csv
```
## References
[1] OpenVINS, "IMU Propagation Derivations," Available: https://docs.openvins.com/propagation.html [Accessed Dec. 30, 2020].\
[2] N. Trawny and S. I. Roumeliotis, "Indirect Kalman Filter for 3D Attitude Estimation," University of Minnesota, Dept. of Comp. Sci. & Eng., Tech. Rep. 2005-002, 2005.\
[3] T. D. Barfoot, *State Estimation for Robotics.* Cambridge University Press, 2017.
