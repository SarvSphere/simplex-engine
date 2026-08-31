cmake_minimum_required(VERSION 3.10)
project(SimplexEngine)

set(CMAKE_CXX_STANDARD 17)

include_directories(include)

add_executable(simplex src/main.cpp src/SimplexSolver.cpp)