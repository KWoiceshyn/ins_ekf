#ifndef INSEKF_UTILITIES_H
#define INSEKF_UTILITIES_H

//#include "ins_ekf_types.h"
#include <eigen3/Eigen/Dense>


namespace ins_ekf{

    Eigen::Vector4d RPYToQuat(const Eigen::Vector3d& rpy);

    Eigen::Matrix3d RPYToRotMat(const Eigen::Vector3d& rpy);

    Eigen::Vector3d QuatToRPY(const Eigen::Vector4d& quat);

    Eigen::Vector4d QuatMultiply(const Eigen::Vector4d& q, const Eigen::Vector4d& p);

    Eigen::Matrix3d SkewSymmetric(const Eigen::Vector3d& vect);

    Eigen::Vector3d RotateVectByQuat(const Eigen::Vector4d& q, const Eigen::Vector3d& v);

    void Print3by3(const Eigen::Matrix3d& mat);

} // ins_ekf

#endif //INSEKF_UTILITIES_H
