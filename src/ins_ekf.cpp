#include "ins_ekf/ins_ekf.h"

#include <iostream>

namespace ins_ekf{

    using MI = MeasurementIndices;
    using SI = StateIndices;

    InsEkf::InsEkf(const Eigen::Vector3d& imu_offset, const Eigen::Vector4d& imu_orientation, const Eigen::Vector3d& gps_offset, Eigen::VectorXd initial_state, double gps_rate):
            imu_offset_{imu_offset},
            imu_orientation_{imu_orientation},
            gps_offset_{gps_offset},
            gps_update_rate_{gps_rate},
            Q_gps_{(Eigen::Matrix<double, 5, 1>() << 0.1, 0.1, 0.1, 0.2, 0.03).finished()},
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

        last_measurements_.setZero();
    }

    void InsEkf::ProcessMeasurement(const std::vector<double>& measurement, double dt) {

        static double t_since_last_gps = 0.0;
        assert(measurement.size() == m_I(MI::NUM_MEASUREMENTS));

        Eigen::Vector3d pqr_meas, axyz_meas;
        pqr_meas << measurement[m_I(MI::p)], measurement[m_I(MI::q)], measurement[m_I(MI::r)];
        axyz_meas << measurement[m_I(MI::ax)], measurement[m_I(MI::ay)], measurement[m_I(MI::az)];

        // transform measurements from sensor frame to vehicle body frame
        pqr_meas = RotateVectByQuat(imu_orientation_, pqr_meas);
        axyz_meas = RotateVectByQuat(imu_orientation_, axyz_meas);

        last_measurements_.segment(m_I(MI::p), 3) = pqr_meas;
        last_measurements_.segment(m_I(MI::ax), 3) = axyz_meas;

        // capture GPS measurements if they are available
        if(t_since_last_gps >= gps_update_rate_){
            const Eigen::VectorXd gps_meas = Eigen::Map<const Eigen::VectorXd>(measurement.data() + 6, 5);
            last_measurements_.segment(m_I(MI::gpsX), 5) = gps_meas;
            t_since_last_gps = 0;
        }else {
            last_measurements_.segment(m_I(MI::gpsX), 5).setConstant(std::nan(""));
            t_since_last_gps += dt;
        }
    }

    void InsEkf::Predict(double dt) { // https://docs.openvins.com/propagation.html

        Eigen::Vector3d w_est = last_measurements_.segment(m_I(MI::p), 3) - bias_g_;
        Eigen::Vector3d a_est = last_measurements_.segment(m_I(MI::ax), 3) - bias_a_;

        // subtract the effect of rotating at offset from vehicle center
        a_est -= SkewSymmetric(w_est) * SkewSymmetric(w_est) * imu_offset_;

        Eigen::Matrix3d w_ss = SkewSymmetric(w_est);
        Eigen::Matrix4d omega = Eigen::Matrix4d::Zero();
        omega.block(0, 0, 3, 3) = -w_ss;
        omega.block(0, 3, 3, 1) = w_est;
        omega.block(3, 0, 1, 3) = -w_est.transpose();

        double angle_rate = w_est.norm();

        // Trawny et al, equation 103, p.12
        Eigen::Matrix4d theta = Eigen::Matrix4d::Identity() + 0.5 * dt * omega;
        if(angle_rate > 1e-10)
            theta = cos(0.5 * angle_rate * dt) * Eigen::Matrix4d::Identity() + (1.0 / angle_rate) * sin(0.5 * angle_rate * dt) * omega;

        // propagate orientation
        orientation_ = theta * orientation_;

        // propagate velocity and position
        Eigen::Matrix3d R_ItoB_hat = RPYToRotMat(QuatToRPY(orientation_)); // rotation matrix estimate from inertial to body frame
        velocity_ = velocity_ - dt * gravity_ + R_ItoB_hat.transpose() * a_est * dt;
        position_ = position_ + dt * velocity_ - 0.5 * gravity_ * dt * dt + 0.5 * R_ItoB_hat.transpose() * a_est * dt * dt;

        // propagate covariance
        double angle = angle_rate * dt;

        // exponential map
        Eigen::Matrix3d R_pert = Eigen::Matrix3d::Identity();
        // SO3 right jacobian
        Eigen::Matrix3d Jr = Eigen::Matrix3d::Identity();

        if(angle > 1e-10){
            Eigen::Vector3d rate_vect = w_est / angle_rate; // normalize

            // Barfoot, equation 6.56, p. 189
            R_pert = cos(angle) * Eigen::Matrix3d::Identity() + (1.0 - cos(angle)) * (rate_vect * rate_vect.transpose()) +
                                     sin(angle) * SkewSymmetric(rate_vect);

            // Barfoot, equation 7.77a, p. 233
            Jr = (sin(angle) / angle) * Eigen::Matrix3d::Identity() + (1.0 - sin(angle) / angle) * (rate_vect * rate_vect.transpose()) -
                 ((1.0 - cos(angle))/angle) * SkewSymmetric(rate_vect);
        }

        Eigen::Matrix<double, s_I(SI::NUM_STATES), s_I(SI::NUM_STATES)> Phi_k = Eigen::Matrix<double, s_I(SI::NUM_STATES), s_I(SI::NUM_STATES)>::Zero();
        Phi_k.block(s_I(SI::roll), s_I(SI::roll), 3, 3) = R_pert;
        Phi_k.block(s_I(SI::x), s_I(SI::roll), 3, 3) = -0.5 * R_ItoB_hat.transpose() * SkewSymmetric(a_est * dt * dt);
        Phi_k.block(s_I(SI::Vx), s_I(SI::roll), 3, 3) = -R_ItoB_hat.transpose() * SkewSymmetric(a_est * dt);

        Phi_k.block(s_I(SI::x), s_I(SI::x), 3, 3).setIdentity();

        Phi_k.block(s_I(SI::x), s_I(SI::Vx), 3, 3) = dt * Eigen::Matrix3d::Identity();
        Phi_k.block(s_I(SI::Vx), s_I(SI::Vx), 3, 3).setIdentity();

        Phi_k.block(s_I(SI::roll), s_I(SI::biasGp), 3, 3) = -R_pert * Jr * dt;
        Eigen::Matrix3d d_w_skew_d_biasGp = SkewSymmetric(Eigen::Vector3d{-1.0, 0.0, 0.0});
        Phi_k.block(s_I(SI::x), s_I(SI::biasGp), 3, 1) = -0.5 * R_ItoB_hat.transpose() * (w_ss*d_w_skew_d_biasGp + d_w_skew_d_biasGp*w_ss) * imu_offset_ * dt * dt;
        Phi_k.block(s_I(SI::Vx), s_I(SI::biasGp), 3, 1) = -R_ItoB_hat.transpose() * (w_ss*d_w_skew_d_biasGp + d_w_skew_d_biasGp*w_ss) * imu_offset_ * dt;
        Eigen::Matrix3d d_w_skew_d_biasGq = SkewSymmetric(Eigen::Vector3d{0.0, -1.0, 0.0});
        Phi_k.block(s_I(SI::x), s_I(SI::biasGq), 3, 1) = -0.5 * R_ItoB_hat.transpose() * (w_ss*d_w_skew_d_biasGq + d_w_skew_d_biasGq*w_ss) * imu_offset_ * dt * dt;
        Phi_k.block(s_I(SI::Vx), s_I(SI::biasGq), 3, 1) = -R_ItoB_hat.transpose() * (w_ss*d_w_skew_d_biasGq + d_w_skew_d_biasGq*w_ss) * imu_offset_ * dt;
        Eigen::Matrix3d d_w_skew_d_biasGr = SkewSymmetric(Eigen::Vector3d{0.0, 0.0, -1.0});
        Phi_k.block(s_I(SI::x), s_I(SI::biasGr), 3, 1) = -0.5 * R_ItoB_hat.transpose() * (w_ss*d_w_skew_d_biasGr + d_w_skew_d_biasGr*w_ss) * imu_offset_ * dt * dt;
        Phi_k.block(s_I(SI::Vx), s_I(SI::biasGr), 3, 1) = -R_ItoB_hat.transpose() * (w_ss*d_w_skew_d_biasGr + d_w_skew_d_biasGr*w_ss) * imu_offset_ * dt;
        Phi_k.block(s_I(SI::biasGp), s_I(SI::biasGp), 3, 3).setIdentity();

        Phi_k.block(s_I(SI::x), s_I(SI::biasAx), 3, 3) = -0.5 * R_ItoB_hat.transpose() * dt * dt;
        Phi_k.block(s_I(SI::Vx), s_I(SI::biasAx), 3, 3) = -R_ItoB_hat.transpose() * dt;
        Phi_k.block(s_I(SI::biasAx), s_I(SI::biasAx), 3, 3).setIdentity();


        Eigen::Matrix<double, s_I(SI::NUM_STATES), 12> Gk = Eigen::Matrix<double, s_I(SI::NUM_STATES), 12>::Zero();
        Gk.block(s_I(SI::roll), m_I(MI::p), 3, 3) = -R_pert * Jr * dt;
        Gk.block(s_I(SI::x), m_I(MI::ax), 3, 3) = -0.5 * R_ItoB_hat.transpose() * dt * dt;
        Gk.block(s_I(SI::Vx), m_I(MI::ax), 3 ,3) = -R_ItoB_hat.transpose() * dt;
        Gk.block(s_I(SI::biasGp), m_I(MI::p) + 6, 3, 3).setIdentity();
        Gk.block(s_I(SI::biasAx), m_I(MI::ax) + 6, 3, 3).setIdentity();

        Eigen::Matrix<double, 12, 12> Qd = Eigen::Matrix<double, 12, 12>::Zero();
        Qd.block(m_I(MI::p), m_I(MI::p), 3, 3) = pow(gyro_noise, 2) * Eigen::Matrix3d::Identity() / dt;
        Qd.block(m_I(MI::ax), m_I(MI::ax), 3, 3) = pow(accel_noise, 2) * Eigen::Matrix3d::Identity() / dt;
        Qd.block(m_I(MI::p) + 6, m_I(MI::p) + 6, 3, 3) = pow(gyro_rw, 2) * Eigen::Matrix3d::Identity() * dt;
        Qd.block(m_I(MI::ax) + 6, m_I(MI::ax) + 6, 3, 3) = pow(accel_rw, 2) * Eigen::Matrix3d::Identity() * dt;

        Eigen::Matrix<double, 15, 15> GQG = Gk * Qd * Gk.transpose();
        GQG = 0.5 * (GQG + GQG.transpose());

        P_ = Phi_k * P_ * Phi_k.transpose() + GQG;

        for(auto i = 0; i < P_.cols(); ++i){
            if(P_(i, i) < 0.0)
                std::cout<<"WARNING: P is not PSD!\n";
        }
    }

    void InsEkf::Update() {

        if(!last_measurements_.hasNaN()){ // need GPS measurements to update

            // measurement model
            std::size_t gps_start = m_I(MI::gpsX);
            Eigen::Matrix<double, 5, 1> h = Eigen::Matrix<double, 5, 1>::Zero();

            Eigen::Vector4d q_BtoI;
            // quaternion of the vehicle's orientation estimate
            q_BtoI << -orientation_[0], -orientation_[1], -orientation_[2], orientation_[3];
            // GPS antenna offset from vehicle center rotated to inertial frame
            Eigen::Vector3d gps_offset_I = RotateVectByQuat(q_BtoI, gps_offset_);
            // expected GPS position measurement based on vehicle position
            h.block(m_I(MI::gpsX) - gps_start, 0, 3, 1) = position_ + gps_offset_I;

            Eigen::Vector3d pqr_est = last_measurements_.segment(m_I(MI::p), 3) - bias_g_;
            // motion of GPS antenna due to rotation at an offset about vehicle center, in body frame
            Eigen::Vector3d gps_rot_offset = SkewSymmetric(pqr_est) * gps_offset_;
            // transform to inertial frame
            Eigen::Vector3d gps_rot_offset_I = RotateVectByQuat(q_BtoI, gps_rot_offset);
            // velocity of GPS antenna in inertial frame
            Eigen::Vector3d gps_velocity_I = velocity_ + gps_rot_offset_I;

            // GPS velocity in inertial XY plane
            double Vxy_norm = gps_velocity_I.segment(0, 2).norm();
            // expeted GPS velocity measurement
            h(m_I(MI::gpsV) - gps_start, 0) = Vxy_norm;
            // expected GPS heading measurement
            h(m_I(MI::gpsPsi) - gps_start, 0) = atan2(gps_velocity_I[1], gps_velocity_I[0]); // atan2(Vy, Vx)

            // measurement model jacobian
            Eigen::Matrix<double, 5, s_I(SI::NUM_STATES)> H = Eigen::Matrix<double, 5, s_I(SI::NUM_STATES)>::Zero();
            H.block(m_I(MI::gpsX) - gps_start, s_I(SI::x), 3, 3) = Eigen::Matrix3d::Identity();

            Eigen::Matrix3d I_3 = Eigen::Matrix3d::Identity();
            // derivative of orientation rotation matrix wrt yaw angle
            Eigen::Matrix3d R_yaw = RPYToRotMat(Eigen::Vector3d{0.0, 0.0, -orientation_[2]});
            // derivative of orientation rotation matrix wrt pitch angle
            Eigen::Matrix3d R_pitch = RPYToRotMat(Eigen::Vector3d{0.0, -orientation_[1], 0.0});
            // derivatives of expected GPS position measurements wrt angle state estimates
            H.block(m_I(MI::gpsX) - gps_start, s_I(SI::roll), 3, 1) = -SkewSymmetric(gps_offset_I) * R_yaw * R_pitch * I_3.col(0); // d gps_offset_I / d roll
            H.block(m_I(MI::gpsX) - gps_start, s_I(SI::pitch), 3, 1) = -SkewSymmetric(gps_offset_I) * R_yaw * I_3.col(1); // d gps_offset_I / d pitch
            H.block(m_I(MI::gpsX) - gps_start, s_I(SI::yaw), 3, 1) = -SkewSymmetric(gps_offset_I) * I_3.col(2); // d gps_offset_I / d yaw

            // derivative of expected GPS velocity measurement wrt velocity state estimates
            H(m_I(MI::gpsV) - gps_start, s_I(SI::Vx)) = gps_velocity_I[0] / Vxy_norm; // Vx / norm(Vxy)
            H(m_I(MI::gpsV) - gps_start, s_I(SI::Vy)) = gps_velocity_I[1] / Vxy_norm; // Vy / norm(Vxy)

            // derivatives of of GPS antenna rotational motion wrt gyro bias state estimates
            Eigen::Vector3d d_gps_roI_d_biasGp = RotateVectByQuat(q_BtoI, SkewSymmetric(Eigen::Vector3d{-1.0, 0.0, 0.0})*gps_offset_);
            Eigen::Vector3d d_gps_roI_d_biasGq = RotateVectByQuat(q_BtoI, SkewSymmetric(Eigen::Vector3d{0.0, -1.0, 0.0})*gps_offset_);
            Eigen::Vector3d d_gps_roI_d_biasGr = RotateVectByQuat(q_BtoI, SkewSymmetric(Eigen::Vector3d{0.0, 0.0, -1.0})*gps_offset_);

            // derivatives of GPS antenna rotational motion wrt angle state estimates
            Eigen::Vector3d d_gps_roI_d_roll = -SkewSymmetric(gps_rot_offset_I) * R_yaw * R_pitch * I_3.col(0); // d gps_rot_offset_I / d roll
            Eigen::Vector3d d_gps_roI_d_pitch = -SkewSymmetric(gps_rot_offset_I) * R_yaw * I_3.col(1); // d gps_rot_offset_I / d pitch
            Eigen::Vector3d d_gps_roI_d_yaw = -SkewSymmetric(gps_rot_offset_I) * I_3.col(2); // d gps_rot_offset_I / d yaw

            // derivative of expected GPS velocity measurement wrt angle state estimates
            H(m_I(MI::gpsV) - gps_start, s_I(SI::roll)) = gps_velocity_I.segment(0, 2).dot(d_gps_roI_d_roll.segment(0, 2)) / Vxy_norm;
            H(m_I(MI::gpsV) - gps_start, s_I(SI::pitch)) = gps_velocity_I.segment(0, 2).dot(d_gps_roI_d_pitch.segment(0, 2)) / Vxy_norm;
            H(m_I(MI::gpsV) - gps_start, s_I(SI::yaw)) = gps_velocity_I.segment(0, 2).dot(d_gps_roI_d_yaw.segment(0, 2)) / Vxy_norm;

            // derivative of expected GPS velocity measurement wrt gyro bias state estimates
            H(m_I(MI::gpsV) - gps_start, s_I(SI::biasGp)) = gps_velocity_I.segment(0, 2).dot(d_gps_roI_d_biasGp.segment(0, 2)) / Vxy_norm;
            H(m_I(MI::gpsV) - gps_start, s_I(SI::biasGq)) = gps_velocity_I.segment(0, 2).dot(d_gps_roI_d_biasGq.segment(0, 2)) / Vxy_norm;
            H(m_I(MI::gpsV) - gps_start, s_I(SI::biasGr)) = gps_velocity_I.segment(0, 2).dot(d_gps_roI_d_biasGr.segment(0, 2)) / Vxy_norm;

            // derivative of expected GPS heading measurement wrt velocity state estimates
            H(m_I(MI::gpsPsi) - gps_start, s_I(SI::Vx)) = -gps_velocity_I[1] / pow(Vxy_norm, 2); // -Vy / norm(Vxy)^2
            H(m_I(MI::gpsPsi) - gps_start, s_I(SI::Vy)) = gps_velocity_I[0] / pow(Vxy_norm, 2); // Vx / norm(Vxy)^2

            // derivative of expected GPS heading measurement wrt angle state estimates
            H(m_I(MI::gpsPsi) - gps_start, s_I(SI::roll)) = (gps_velocity_I[0] * d_gps_roI_d_roll[1] - gps_velocity_I[1] * d_gps_roI_d_roll[0]) / pow(Vxy_norm, 2); // (Vx*dVy/droll - Vy*dVx/droll) / norm^2
            H(m_I(MI::gpsPsi) - gps_start, s_I(SI::pitch)) = (gps_velocity_I[0] * d_gps_roI_d_pitch[1] - gps_velocity_I[1] * d_gps_roI_d_pitch[0]) / pow(Vxy_norm, 2);
            H(m_I(MI::gpsPsi) - gps_start, s_I(SI::yaw)) = (gps_velocity_I[0] * d_gps_roI_d_yaw[1] - gps_velocity_I[1] * d_gps_roI_d_yaw[0]) / pow(Vxy_norm, 2);

            // derivative of expected GPS heading measurement wrt gyro bias state estimates
            H(m_I(MI::gpsPsi) - gps_start, s_I(SI::biasGp)) = (gps_velocity_I[0] * d_gps_roI_d_biasGp[1] - gps_velocity_I[1] * d_gps_roI_d_biasGp[0]) / pow(Vxy_norm, 2);
            H(m_I(MI::gpsPsi) - gps_start, s_I(SI::biasGq)) = (gps_velocity_I[0] * d_gps_roI_d_biasGq[1] - gps_velocity_I[1] * d_gps_roI_d_biasGq[0]) / pow(Vxy_norm, 2);
            H(m_I(MI::gpsPsi) - gps_start, s_I(SI::biasGr)) = (gps_velocity_I[0] * d_gps_roI_d_biasGr[1] - gps_velocity_I[1] * d_gps_roI_d_biasGr[0]) / pow(Vxy_norm, 2);

            // EKF update
            const Eigen::Matrix<double, 5, 5> Q_gps_diag = Q_gps_.asDiagonal();
            // covariance of residual
            Eigen::Matrix<double, 5, 5> S = H * P_ * H.transpose() + Q_gps_diag;

            Eigen::Matrix<double, 5, 5> S_inverse = Eigen::Matrix<double, 5, 5>::Identity();
            S.selfadjointView<Eigen::Upper>().llt().solveInPlace(S_inverse);
            // Kalman gain = P * H'/ S
            Eigen::Matrix<double, s_I(SI::NUM_STATES), 5> K = P_ * H.transpose() * S_inverse;

            Eigen::Matrix<double, 5, 1> residual = last_measurements_.segment(m_I(MI::gpsX), 5) - h;

            // update the states
            Eigen::Matrix<double, s_I(SI::NUM_STATES), 1> dx = K * residual;
            Eigen::Vector4d dq {0.5 * dx[s_I(SI::roll)], 0.5 * dx[s_I(SI::pitch)], 0.5 * dx[s_I(SI::yaw)], 1.0}; // small angle approximation
            dq /= dq.norm();
            orientation_ = QuatMultiply(dq, orientation_);
            orientation_ /= orientation_.norm();
            position_ += dx.block(s_I(SI::x), 0, 3, 1);
            velocity_ += dx.block(s_I(SI::Vx), 0, 3, 1);
            bias_g_ += dx.block(s_I(SI::biasGp), 0, 3, 1);
            bias_a_ += dx.block(s_I(SI::biasAx), 0, 3, 1);

            // update the covariance
            Eigen::Matrix<double, s_I(SI::NUM_STATES), s_I(SI::NUM_STATES)> I_KH =
                    Eigen::Matrix<double, s_I(SI::NUM_STATES), s_I(SI::NUM_STATES)>::Identity() - K * H;

            P_ = I_KH * P_ * I_KH.transpose() + K * Q_gps_diag * K.transpose();
            P_ = 0.5 * (P_ + P_.transpose());
            //P_ -= K * H * P_;

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

