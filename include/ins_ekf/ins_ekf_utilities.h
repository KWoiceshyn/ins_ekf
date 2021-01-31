#ifndef INSEKF_UTILITIES_H
#define INSEKF_UTILITIES_H

#include <eigen3/Eigen/Dense>

namespace ins_ekf{

    // Euler roll,pitch,yaw angles to quaternion
    // https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
    Eigen::Vector4d RPYToQuat(const Eigen::Vector3d& rpy);

    // Euler roll,pitch,yaw angles to rotation matrix
    Eigen::Matrix3d RPYToRotMat(const Eigen::Vector3d& rpy);

    // Quaternion to Euler roll,pitch,yaw angles
    // https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
    Eigen::Vector3d QuatToRPY(const Eigen::Vector4d& quat);

    // Multiply two quaternions
    Eigen::Vector4d QuatMultiply(const Eigen::Vector4d& q, const Eigen::Vector4d& p);

    // Return a 3x3 skew symmetric matrix from an R3 vector
    Eigen::Matrix3d SkewSymmetric(const Eigen::Vector3d& vect);

    // Rotate an R3 vector by a quaternion
    Eigen::Vector3d RotateVectByQuat(const Eigen::Vector4d& q, const Eigen::Vector3d& v);

    // Print the contents of a 3x3 matrix
    void Print3by3(const Eigen::Matrix3d& mat);

} // namespace ins_ekf

#endif //INSEKF_UTILITIES_H
