#include "ins_ekf/ins_ekf.h"

#include <fstream>
#include <sstream>
#include <iostream>

using namespace ins_ekf;

int main(int argc, char** argv){

    if(argc < 3){
        std::cerr << "Usage: " << argv[0] << "--mf measurements_file \n";
        return 1;
    }

    Eigen::Vector3d imu_offset(0.5, 0.5, 0.5);
    Eigen::Vector4d imu_orientation(0.0, 0.0, 0.0, 1.0);
    Eigen::Vector3d gps_offset(1.0, 1.0, 1.0);

    Eigen::Matrix<double, s_I(StateIndices::NUM_STATES), 1> initial_state;

    initial_state << 0.0, 0.0, 0.0, // RPY
            0.0, 0.0, 0.0, // xyz position
            0.0, 0.0, 0.0, // Vxyz
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0;

    ins_ekf::InsEkf test_ekf(imu_offset, imu_orientation, gps_offset, initial_state, 0.1);

    double dt = 0.01;

    std::ifstream meas_file;
    std::ofstream outstates_file;
    std::ofstream outcov_file;
    meas_file.open(std::string(argv[2]));
    outstates_file.open("../states_out.csv");
    outcov_file.open("../cov_out.csv");
    if(!meas_file.is_open()){
        std::cerr << "Measurements file not opened.\n";
        return 1;
    }

    double time = 0.0;

    std::string line;
    while(getline(meas_file, line)){
        std::istringstream iss(line);
        std::vector<double> meas;
        std::string line_stream;
        std::string::size_type sz;
        while(getline(iss, line_stream, ',')){
            meas.push_back(stod(line_stream, &sz));
        }

        test_ekf.ProcessMeasurement(meas, dt);
        test_ekf.Predict(dt);
        test_ekf.Update();

        auto states = test_ekf.GetStates();
        auto cov = test_ekf.GetCovariance();

        time += dt;

        outstates_file << time << ",";
        outcov_file << time << ",";
        for(auto i = 0; i < states.size(); ++i){
            outstates_file << states[i];
            outcov_file << cov[i];
            if(i < states.size() - 1){
                outstates_file << ",";
                outcov_file << ",";
            }
        }
        outstates_file << "\n";
        outcov_file << "\n";

    }

    meas_file.close();
    outstates_file.close();
    outcov_file.close();

    return 0;
}

