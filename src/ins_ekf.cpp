#include <iostream>

#include "../include/ins_ekf/ins_ekf.h"

namespace ins_ekf{

    using MI = MeasurementIndices;
    using SI = StateIndices;

    InsEkf::InsEkf(const Eigen::Vector3d& imu_offset, const Eigen::Vector4d& imu_orientation, const Eigen::Vector3d& gps_offset, Eigen::VectorXd initial_state):
            imu_offset_{imu_offset},
            imu_orientation_{imu_orientation},
            gps_offset_{gps_offset},
            gravity_{(Eigen::Matrix<double, 3, 1>() << 0.0, 0.0, 9.81).finished()},
            is_stationary_{false}
    {

        assert(initial_state.size() == s_I(SI::NUM_STATES));
        orientation_ = RPYToQuat(initial_state.segment(s_I(SI::roll), 3));
        position_ = initial_state.segment(s_I(SI::x), 3);
        velocity_ = initial_state.segment(s_I(SI::Vx), 3);
        bias_g_.setZero();
        bias_a_.setZero();

        // initialize state covariance matrix
        P_.setZero();
        P_.block(s_I(SI::x), s_I(SI::x), 3, 3) = 10.0 * Eigen::Matrix3d::Identity(); // xyz
        P_.block(s_I(SI::roll), s_I(SI::roll), 2 , 2) = 8e-3 * Eigen::Matrix2d::Identity(); // roll, pitch
        P_(s_I(SI::yaw), s_I(SI::yaw)) = 3e-2; // yaw
        P_.block(s_I(SI::Vx), s_I(SI::Vx), 3 , 3) = 0.1 * Eigen::Matrix3d::Identity(); // Vxyz
        P_.block(s_I(SI::biasGp), s_I(SI::biasGp), 3, 3) = 1e-7 * Eigen::Matrix3d::Identity(); // biasGpqr
        P_.block(s_I(SI::biasAx), s_I(SI::biasAx), 3, 3) = 1e-8 * Eigen::Matrix3d::Identity(); // biasAxyz

        // GPS noise values
        Q_gps_ << 0.1, 0.1, 0.1, 0.2, 0.03;

        double accel_noise = 1e-3; // m/s^2/hz^0.5
        double accel_rw = 1e-3; // m/s^2/hz^0.5
        double gyro_noise = 1.3e-4; // rad/s/hz^0.5
        double gyro_rw = 1e-4; // rad/s/hz^0.5

        Q_imu_.setZero();
        Q_imu_.block(0, 0, 3, 3) = pow(gyro_noise, 2) * Eigen::Matrix3d::Identity();
        Q_imu_.block(3, 3, 3, 3) = pow(accel_noise, 2) * Eigen::Matrix3d::Identity();
        Q_imu_.block(6, 6, 3, 3) = pow(gyro_rw, 2) * Eigen::Matrix3d::Identity();
        Q_imu_.block(9, 9, 3, 3) = pow(accel_rw, 2) * Eigen::Matrix3d::Identity();

        last_measurements_.setZero();

    }

    void InsEkf::ProcessMeasurement(const std::vector<double>& measurement, double dt) {

        static double t_since_last_gps = 0.0;
        assert(measurement.size() == m_I(MI::NUM_MEASUREMENTS));

        const Eigen::VectorXd imu_meas = Eigen::Map<const Eigen::VectorXd>(measurement.data(), 6);
        last_measurements_.segment(m_I(MI::p), 6) = imu_meas;

        if(t_since_last_gps >= 0.1){
            const Eigen::VectorXd gps_meas = Eigen::Map<const Eigen::VectorXd>(measurement.data() + 6, 5);
            last_measurements_.segment(m_I(MI::gpsX), 5) = gps_meas;
            t_since_last_gps = 0;
        }else {
            last_measurements_.segment(m_I(MI::gpsX), 5).setConstant(std::nan(""));
            t_since_last_gps += dt;
        }
    }

    void InsEkf::Predict(double dt) { // https://docs.openvins.com/propagation.html

        // propagate states
        Eigen::Vector3d w_est = last_measurements_.segment(m_I(MI::p), 3) - bias_g_;
        Eigen::Vector3d a_est = last_measurements_.segment(m_I(MI::ax), 3) - bias_a_;
        Eigen::Matrix3d w_ss = SkewSymmetric(w_est);
        Eigen::Matrix4d omega = Eigen::Matrix4d::Zero();
        omega.block(0, 0, 3, 3) = -w_ss;
        omega.block(0, 3, 3, 1) = w_est;
        omega.block(3, 0, 1, 3) = -w_est.transpose();

        double angle_rate = w_est.norm();

        Eigen::Matrix4d theta = Eigen::Matrix4d::Identity() + 0.5 * dt * omega;
        if(angle_rate > 1e-10)
            theta = cos(0.5 * angle_rate * dt) * Eigen::Matrix4d::Identity() + (1.0 / angle_rate) * sin(0.5 * angle_rate * dt) * omega;

        orientation_ = theta * orientation_;

        Eigen::Matrix3d R_ItoB_hat = RPYToRotMat(QuatToRPY(orientation_));
        velocity_ = velocity_ - dt * gravity_ + R_ItoB_hat.transpose() * (last_measurements_.segment(m_I(MI::ax), 3) - bias_a_) * dt;
        position_ = position_ + dt * velocity_ - 0.5 * gravity_ * dt * dt + 0.5 * R_ItoB_hat.transpose() * (last_measurements_.segment(m_I(MI::ax), 3) - bias_a_) * dt * dt;

        // propagate covariance
        double angle = angle_rate * dt; // TODO : this may need to be negative?

        // exponential map
        Eigen::Matrix3d R_pert = Eigen::Matrix3d::Identity();
        // right jacobian
        Eigen::Matrix3d Jr = Eigen::Matrix3d::Identity();

        if(angle > 1e-10){
            Eigen::Vector3d rate_vect = w_est / angle_rate;

            R_pert = cos(angle) * Eigen::Matrix3d::Identity() + (1.0 - cos(angle)) * rate_vect * rate_vect.transpose() +
                                     sin(angle) * SkewSymmetric(rate_vect);

            Jr = (sin(angle) / angle) * Eigen::Matrix3d::Identity() + (1.0 - sin(angle) / angle) * rate_vect * rate_vect.transpose() -
                 ((1.0 - cos(angle))/angle) * SkewSymmetric(rate_vect);
        }

        Eigen::Matrix<double, 15, 15> Phi_k = Eigen::Matrix<double, 15, 15>::Zero();
        Phi_k.block(0, 0, 3, 3) = R_pert;
        Phi_k.block(3, 0, 3, 3) = -0.5 * R_ItoB_hat.transpose() * SkewSymmetric(a_est * dt * dt);
        Phi_k.block(6, 0, 3, 3) = -R_ItoB_hat.transpose() * SkewSymmetric(a_est * dt);

        Phi_k.block(3, 3, 3, 3).setIdentity();

        Phi_k.block(3, 6, 3, 3) = dt * Eigen::Matrix3d::Identity();
        Phi_k.block(6, 6, 3, 3).setIdentity();

        Phi_k.block(0, 9, 3, 3) = -R_pert * Jr * dt;
        Phi_k.block(9, 9, 3, 3).setIdentity();

        Phi_k.block(12, 3, 3, 3) = -0.5 * R_ItoB_hat.transpose() * dt * dt;
        Phi_k.block(12, 6, 3, 3) = -R_ItoB_hat.transpose() * dt;
        Phi_k.block(12, 12, 3, 3).setIdentity();


        Eigen::Matrix<double, 15, 12> Gk = Eigen::Matrix<double, 15, 12>::Zero();
        Gk.block(0, 0, 3, 3) = -R_pert * Jr * dt;
        Gk.block(3, 3, 3, 3) = -0.5 * R_ItoB_hat.transpose() * dt * dt;
        Gk.block(6, 3, 3 ,3) = -R_ItoB_hat.transpose() * dt;
        Gk.block(9, 6, 3, 3).setIdentity();
        Gk.block(12, 9, 3, 3).setIdentity();

        P_ = Phi_k * P_ * Phi_k.transpose() + Gk * Q_imu_ * Gk.transpose();
    }

    std::vector<double> InsEkf::GetStates() const {

        std::vector<double> states;
        Eigen::Vector3d orientation = QuatToRPY(orientation_);

        for(auto i = 0; i < 3; ++i)
            states.push_back(position_[i]);
        for(auto i = 0; i < 3; ++i)
            states.push_back(orientation[i]);
        for(auto i = 0; i < 3; ++i)
            states.push_back(velocity_[i]);
        for(auto i = 0; i < 3; ++i)
            states.push_back(bias_g_[i]);
        for(auto i = 0; i < 3; ++i)
            states.push_back(bias_a_[i]);

        return states;
    }

    std::vector<double> InsEkf::GetCovariance() const {

        std::vector<double> cov;

        for(auto i = 0; i < s_I(SI::NUM_STATES); ++i)
            cov.push_back(P_(i, i)); // diagonals

        return cov;
    }

}

