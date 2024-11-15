# Copyright 2024 Thespica
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.0)
project(play_mpi C)

set(CMAKE_C_STANDARD 11)

set(MPI_EXECUTABLE_SUFFIX ".mpich")
find_package(MPI REQUIRED)

add_executable(play_mpi main.c)
target_link_libraries(play_mpi PRIVATE ${MPI_C_LIBRARIES}:)
