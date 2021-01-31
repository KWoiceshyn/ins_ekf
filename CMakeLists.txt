cmake_minimum_required (VERSION 3.11)
set(CMAKE_CXX_STANDARD 14)
project(ins_ekf)

add_library(ins_ekf src/ins_ekf.cpp)
target_include_directories(ins_ekf PUBLIC include)

add_executable(main src/main.cpp)
target_link_libraries(main ins_ekf)
