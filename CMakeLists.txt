cmake_minimum_required (VERSION 3.11)
set(CMAKE_CXX_STANDARD 14)
project(ins_ekf)

set(SOURCES
        src/ins_ekf_utilities.cpp
        src/ins_ekf.cpp
        )

add_library(ins_ekf ${SOURCES})
target_include_directories(ins_ekf PUBLIC include)

add_executable(main src/main.cpp)
target_link_libraries(main ins_ekf)
