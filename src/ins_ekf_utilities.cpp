#include "../include/ins_ekf/ins_ekf_utilities.h"

namespace ins_ekf{

    Eigen::Vector4d RPYToQuat(const Eigen::Vector3d& rpy){ // https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles

        Eigen::Vector4d quat;
        double cr_2   = cos(rpy[0] / 2);
        double sr_2   = sin(rpy[0] / 2);
        double cp_2   = cos(rpy[1] / 2);
        double sp_2   = sin(rpy[1] / 2);
        double cy_2   = cos(rpy[2] / 2);
        double sy_2   = sin(rpy[2] / 2);

        quat[0] = sr_2 * cp_2 * cy_2 - cr_2 * sp_2 * sy_2; // x
        quat[1] = cr_2 * sp_2 * cy_2 + sr_2 * cp_2 * sy_2; // y
        quat[2] = cr_2 * cp_2 * sy_2 - sr_2 * sp_2 * cy_2; // z
        quat[3] = cr_2 * cp_2 * cy_2 + sr_2 * sp_2 * sy_2; // w

        quat /= quat.norm();

        return quat;
    }

    Eigen::Matrix3d RPYToRotMat(const Eigen::Vector3d& rpy){

        Eigen::Matrix3d rot_I_to_B; // inertial frame to body frame
        double cr   = cos(rpy[0]);
        double sr   = sin(rpy[0]);
        double cp   = cos(rpy[1]);
        double sp   = sin(rpy[1]);
        double cy   = cos(rpy[2]);
        double sy   = sin(rpy[2]);

        rot_I_to_B << cp*cy, cp*sy, -sp,
                sr*sp*cy - cr*sy, sr*sp*sy + cr*cy, sr*cp,
                cr*sp*cy + sr*sy, cr*sp*sy - sr*cy, cr*cp;

        return rot_I_to_B;
    }

    Eigen::Vector3d QuatToRPY(const Eigen::Vector4d& quat){ // https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles

        // roll (body x axis)
        double sr_cp = 2.0 * (quat[3] * quat[0] + quat[1] * quat[2]);
        double cr_cp = 1.0 - 2.0 * (quat[0] * quat[0] + quat[1] * quat[1]);
        double roll = atan2(sr_cp, cr_cp);

        // pitch (y axis of pre-roll frame)
        double sp = 2.0 * (quat[3] * quat[1] - quat[2] * quat[0]);
        double pitch = 0.0;
        if (std::fabs(sp) > 0.98) // unlikely to be vertical pitch
            pitch = sp > 0 ? M_PI_2 : -M_PI_2;
        else
            pitch = asin(sp);

        // yaw (inertial frame z axis)
        double sy_cp = 2.0 * (quat[3] * quat[2] + quat[0] * quat[1]);
        double cy_cp = 1.0 - 2.0 * (quat[1] * quat[1] + quat[2] * quat[2]);
        double yaw = atan2(sy_cp, cy_cp);

        Eigen::Vector3d rpy {roll, pitch, yaw};

        return rpy;
    }

    Eigen::Vector4d QuatMultiply(const Eigen::Vector4d& q, const Eigen::Vector4d& p){

        Eigen::Vector4d quat;

        quat[0] = q[3]*p[0] + q[2]*p[1] - q[1]*p[2] + q[0]*p[3];
        quat[1] = -q[2]*p[0] + q[3]*p[1] + q[0]*p[2] + q[1]*p[3];
        quat[2] = q[1]*p[0] - q[0]*p[1] + q[3]*p[2] + q[2]*p[3];
        quat[3] = -q[0]*p[0] - q[1]*p[1] - q[2]*p[2] + q[3]*p[3];

        return quat;
    }

    Eigen::Matrix3d SkewSymmetric(const Eigen::Vector3d& vect){

        Eigen::Matrix3d ss;
        ss << 0, -vect[2], vect[1],
                vect[2], 0, -vect[0],
                -vect[1], vect[0], 0;
        return ss;
    }

    Eigen::Vector3d RotateVectByQuat(const Eigen::Vector4d& q, const Eigen::Vector3d& v){

        if(v.norm() < 1e-10)
            return v;

        Eigen::Vector4d q_inv, v_quat;
        q_inv << -q[0], -q[1], -q[2], q[3];
        v_quat << v[0], v[1], v[2], 0.0;

        Eigen::Vector4d result = QuatMultiply(v_quat, q_inv);
        result = QuatMultiply(q, result);
        return result.segment(0, 3);
    }

    void Print3by3(const Eigen::Matrix3d& mat){

        printf("%4.4lf\t\t%4.4lf\t\t%4.4lf\n", mat(0,0), mat(0,1), mat(0,2));
        printf("%4.4lf\t\t%4.4lf\t\t%4.4lf\n", mat(1,0), mat(1,1), mat(1,2));
        printf("%4.4lf\t\t%4.4lf\t\t%4.4lf\n", mat(2,0), mat(2,1), mat(2,2));

    }

} // ins_ekf