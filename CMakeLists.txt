cmake_minimum_required(VERSION 3.10)
project(matrix_calculator)

set(CMAKE_CXX_STANDARD 17)
include_directories(head)

add_executable(matrix_calculator main.cpp source/Rational.cpp)

