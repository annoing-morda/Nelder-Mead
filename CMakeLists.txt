cmake_minimum_required(VERSION 2.8.4)
project(neldermead)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")

set(SOURCE_FILES main.cpp neldermead.cpp)
add_executable(neldermead ${SOURCE_FILES})