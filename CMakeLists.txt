cmake_minimum_required(VERSION 2.6)
cmake_policy(SET CMP0015 NEW)
project(cpdeflec)

set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -std=c99 -O2")

#include_directories()
#link_directories()

add_executable(cpdeflec
	cpdeflec.c
	commath.c
	camera.c
	pattern.c
	fitting.c
	)

target_link_libraries(cpdeflec gsl gslcblas tiff m)
