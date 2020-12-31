#ifndef INSEKF_UTILITIES_H
#define INSEKF_UTILITIES_H

//#include "ins_ekf_types.h"
#include <eigen3/Eigen/Dense>


namespace ins_ekf{

    Eigen::Vector4d RPYToQuat(const Eigen::Vector3d& rpy);

    Eigen::Matrix3d RPYToRotMat(const Eigen::Vector3d& rpy);

    Eigen::Vector3d QuatToRPY(const Eigen::Vector4d& quat);

    Eigen::Matrix3d SkewSymmetric(const Eigen::Vector3d& vect);

} // ins_ekf

#endif //INSEKF_UTILITIES_H
