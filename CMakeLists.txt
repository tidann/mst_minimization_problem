cmake_minimum_required(VERSION 3.10)
project(MST_Project)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

add_executable(mst_solver mst_solver.cpp)
