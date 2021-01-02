#ifndef INSEKF_INSEKFTYPES_H
#define INSEKF_INSEKFTYPES_H

#include <vector>

namespace ins_ekf {

    enum class MeasurementIndices {
        p, // IMU angular rates
        q,
        r,
        ax, // IMU accelerations
        ay,
        az,
        gpsX, // GPS position
        gpsY,
        gpsZ,
        gpsV, // GPS velocity over ground
        gpsPsi, // GPS track over ground (heading)
        NUM_MEASUREMENTS
    };

    enum class StateIndices {
        roll, // attitude
        pitch,
        yaw,
        x, // inertial frame position
        y,
        z,
        Vx, // inertial frame velocities
        Vy,
        Vz,
        biasGp, // gyro bias
        biasGq,
        biasGr,
        biasAx, // accelerometer bias
        biasAy,
        biasAz,
        NUM_STATES
    };

    constexpr std::size_t m_I(MeasurementIndices m){return static_cast<std::size_t>(m);};
    constexpr std::size_t s_I(StateIndices s){return static_cast<std::size_t>(s);};

} // ins_ekf

#endif //INSEKF_INSEKFTYPES_H