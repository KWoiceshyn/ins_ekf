#ifndef INSEKF_INSEKF_H
#define INSEKF_INSEKF_H

#include "ins_ekf/ins_ekf_types.h"
#include "ins_ekf/ins_ekf_utilities.h"

#include <vector>

namespace ins_ekf{

class InsEkf{

public:
    InsEkf(const Eigen::Vector3d& imu_offset, const Eigen::Vector4d& imu_orientation, const Eigen::Vector3d& gps_offset,
           Eigen::VectorXd initial_state, double gps_rate);

    void ProcessMeasurement(const std::vector<double>& measurement, double dt);

    void Predict(double dt);

    void Update();

    std::vector<double> GetStates() const;
    std::vector<double> GetCovariance() const;

private:

    // state estimates
    Eigen::Matrix<double, 4, 1> orientation_;
    Eigen::Matrix<double, 3, 1> position_;
    Eigen::Matrix<double, 3, 1> velocity_;
    Eigen::Matrix<double, 3, 1> bias_g_;
    Eigen::Matrix<double, 3, 1> bias_a_;


    Eigen::Matrix<double, s_I(StateIndices ::NUM_STATES), s_I(StateIndices::NUM_STATES)> P_; // state covariance

    const double accel_noise = 1e-3; // m/s^2/hz^0.5
    const double accel_rw = 1e-3; // m/s^2/hz^0.5
    const double gyro_noise = 1.3e-4; // rad/s/hz^0.5
    const double gyro_rw = 1e-4; // rad/s/hz^0.5
    const double gps_update_rate_; // sec

    const Eigen::Matrix<double, 5, 1> Q_gps_; // GPS measurement noises

    Eigen::Matrix<double, m_I(MeasurementIndices::NUM_MEASUREMENTS), 1> last_measurements_;

    const Eigen::Matrix<double, 3, 1> gravity_; // {0, 0, 9.81}

    bool is_stationary_;

    const Eigen::Vector3d imu_offset_; // xyz offset of IMU from vehicle rotational center
    const Eigen::Vector4d imu_orientation_; // orientation of IMU in vehicle body frame
    const Eigen::Vector3d gps_offset_; // xyz offset of GPS antenna from vehicle rotational center
};

} // ins_ekf


#endif //INSEKF_INSEKF_H
