cmake_minimum_required(VERSION 3.0)
project(play_mpi C)

set(CMAKE_C_STANDARD 11)

set(MPI_EXECUTABLE_SUFFIX ".mpich")
find_package(MPI REQUIRED)

add_executable(play_mpi main.c)
target_link_libraries(play_mpi PRIVATE ${MPI_C_LIBRARIES}:)
