cmake_minimum_required(VERSION 3.7)
project(MPI_Test_Programs)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_C_STANDARD 11)

set(SOURCE_FILES main.c)
add_executable(MPI_Test_Programs ${SOURCE_FILES})