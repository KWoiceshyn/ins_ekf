cmake_minimum_required (VERSION 3.11)
set(CMAKE_CXX_STANDARD 14)
project(INSEkf)

set(SOURCES
		"src/ins_ekf_utilities.cpp"
		"src/ins_ekf.cpp"
		"src/main.cpp"

		)
set(HEADERS
		"include/ins_ekf/ins_ekf_types.h"
		"include/ins_ekf/ins_ekf_utilities.h"
		"include/ins_ekf/ins_ekf.h"
		)

add_executable(${PROJECT_NAME}
		${SOURCES}
		${HEADERS}
		)
