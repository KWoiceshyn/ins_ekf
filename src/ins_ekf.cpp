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
        Eigen::Matrix<double, 5, 1> gps_noise;
        gps_noise << 0.1, 0.1, 0.1, 0.2, 0.03;
        Q_gps_ = gps_noise.asDiagonal();

        double accel_noise = 1e-3; // m/s^2/hz^0.5
        double accel_rw = 1e-3; // m/s^2/hz^0.5
        double gyro_noise = 1.3e-4; // rad/s/hz^0.5
        double gyro_rw = 1e-4; // rad/s/hz^0.5

        Q_imu_.setZero();
        Q_imu_.block(m_I(MI::p), m_I(MI::p), 3, 3) = pow(gyro_noise, 2) * Eigen::Matrix3d::Identity();
        Q_imu_.block(m_I(MI::ax), m_I(MI::ax), 3, 3) = pow(accel_noise, 2) * Eigen::Matrix3d::Identity();
        Q_imu_.block(m_I(MI::p) + 6, m_I(MI::p) + 6, 3, 3) = pow(gyro_rw, 2) * Eigen::Matrix3d::Identity();
        Q_imu_.block(m_I(MI::ax) + 6, m_I(MI::ax) + 6, 3, 3) = pow(accel_rw, 2) * Eigen::Matrix3d::Identity();

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
        Phi_k.block(s_I(SI::roll), s_I(SI::roll), 3, 3) = R_pert;
        Phi_k.block(s_I(SI::x), s_I(SI::roll), 3, 3) = -0.5 * R_ItoB_hat.transpose() * SkewSymmetric(a_est * dt * dt);
        Phi_k.block(s_I(SI::Vx), s_I(SI::roll), 3, 3) = -R_ItoB_hat.transpose() * SkewSymmetric(a_est * dt);

        Phi_k.block(s_I(SI::x), s_I(SI::x), 3, 3).setIdentity();

        Phi_k.block(s_I(SI::x), s_I(SI::Vx), 3, 3) = dt * Eigen::Matrix3d::Identity();
        Phi_k.block(s_I(SI::Vx), s_I(SI::Vx), 3, 3).setIdentity();

        Phi_k.block(s_I(SI::roll), s_I(SI::biasGp), 3, 3) = -R_pert * Jr * dt;
        Phi_k.block(s_I(SI::biasGp), s_I(SI::biasGp), 3, 3).setIdentity();

        Phi_k.block(s_I(SI::x), s_I(SI::biasAx), 3, 3) = -0.5 * R_ItoB_hat.transpose() * dt * dt;
        Phi_k.block(s_I(SI::Vx), s_I(SI::biasAx), 3, 3) = -R_ItoB_hat.transpose() * dt;
        Phi_k.block(s_I(SI::biasAx), s_I(SI::biasAx), 3, 3).setIdentity();


        Eigen::Matrix<double, 15, 12> Gk = Eigen::Matrix<double, 15, 12>::Zero();
        Gk.block(s_I(SI::roll), m_I(MI::p), 3, 3) = -R_pert * Jr * dt;
        Gk.block(s_I(SI::x), m_I(MI::ax), 3, 3) = -0.5 * R_ItoB_hat.transpose() * dt * dt;
        Gk.block(s_I(SI::Vx), m_I(MI::ax), 3 ,3) = -R_ItoB_hat.transpose() * dt;
        Gk.block(s_I(SI::biasGp), m_I(MI::p) + 6, 3, 3).setIdentity();
        Gk.block(s_I(SI::biasAx), m_I(MI::ax) + 6, 3, 3).setIdentity();

        Eigen::Matrix<double, 12, 12> Qd = Q_imu_;
        Qd.block(m_I(MI::p), m_I(MI::p), 6, 6) *= (1.0 / dt);
        Qd.block(m_I(MI::p) + 6, m_I(MI::p) + 6, 6, 6) *= dt;

        P_ = Phi_k * P_ * Phi_k.transpose() + Gk * Qd * Gk.transpose();
    }

    void InsEkf::Update() {

        if(!last_measurements_.hasNaN()){ // need GPS measurements to update

            // measurement model
            std::size_t gps_start = m_I(MI::gpsX);
            Eigen::Matrix<double, 5, 1> h = Eigen::Matrix<double, 5, 1>::Zero();
            h.block(m_I(MI::gpsX) - gps_start, 0, 3, 1) = position_; // TODO : incorporate GPS sensor offset
            h(m_I(MI::gpsV) - gps_start, 0) = sqrt(pow(velocity_[0], 2) + pow(velocity_[1], 2)); // norm(Vxy)
            h(m_I(MI::gpsPsi) - gps_start, 0) = atan2(velocity_[1], velocity_[0]);

            // measurement model jacobian
            Eigen::Matrix<double, 5, s_I(SI::NUM_STATES)> H = Eigen::Matrix<double, 5, s_I(SI::NUM_STATES)>::Zero();
            H.block(0, s_I(SI::x), 3, 3) = Eigen::Matrix3d::Identity(); // TODO : incorporate GPS sensor offset
            H(m_I(MI::gpsV) - gps_start, s_I(SI::Vx)) = velocity_[0] / h(m_I(MI::gpsV) - gps_start, 0); // Vx / norm(Vxy)
            H(m_I(MI::gpsV) - gps_start, s_I(SI::Vy)) = velocity_[1] / h(m_I(MI::gpsV) - gps_start, 0); // Vy / norm(Vxy)
            H(m_I(MI::gpsPsi) - gps_start, s_I(SI::Vx)) = -velocity_[1] / pow(h(m_I(MI::gpsV) - gps_start, 0), 2); // -Vy / norm(Vxy)^2
            H(m_I(MI::gpsPsi) - gps_start, s_I(SI::Vy)) = velocity_[0] / pow(h(m_I(MI::gpsV) - gps_start, 0), 2); // Vx / norm(Vxy)^2

            Eigen::Matrix<double, 5, 5> S = H * P_ * H.transpose() + Q_gps_;

            // Kalman gain = P * H'/ S
            Eigen::Matrix<double, 5, 5> S_inverse;
            S.selfadjointView<Eigen::Upper>().llt().solveInPlace(S_inverse);
            Eigen::Matrix<double, s_I(SI::NUM_STATES), 5> K = P_ * H.transpose() * S_inverse;

            Eigen::Matrix<double, 5, 1> residual = last_measurements_.segment(m_I(MI::gpsX), 5) - h;

            // update the states
            Eigen::Matrix<double, s_I(SI::NUM_STATES), 1> dx = K * residual;
            Eigen::Vector4d dq {0.5 * dx[s_I(SI::roll)], 0.5 * dx[s_I(SI::pitch)], 0.5 * dx[s_I(SI::yaw)], 1.0}; // small angle approximation
            dq /= dq.norm();
            orientation_ = QuatMultiply(dq, orientation_);
            position_ += dx.block(s_I(SI::x), 0, 3, 1);
            velocity_ += dx.block(s_I(SI::Vx), 0, 3, 1);
            bias_g_ += dx.block(s_I(SI::biasGp), 0, 3, 1);
            bias_a_ += dx.block(s_I(SI::biasAx), 0, 3, 1);

            // update the covariance
            Eigen::Matrix<double, s_I(SI::NUM_STATES), s_I(SI::NUM_STATES)> I_KH =
                    Eigen::Matrix<double, s_I(SI::NUM_STATES), s_I(SI::NUM_STATES)>::Identity() - K * H;

            P_ = I_KH * P_ * I_KH.transpose() + K * Q_gps_ * K.transpose();

        }
    }

    std::vector<double> InsEkf::GetStates() const {

        std::vector<double> states;
        Eigen::Vector3d orientation = QuatToRPY(orientation_);

        for(auto i = 0; i < 3; ++i)
            states.push_back(orientation[i]);
        for(auto i = 0; i < 3; ++i)
            states.push_back(position_[i]);
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

