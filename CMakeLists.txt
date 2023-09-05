cmake_minimum_required(VERSION 3.16)
project(kurs3chm)

set(CMAKE_CXX_STANDARD 14)
#Set Linker flags
find_library(SOMELIB libstdc++)
find_library(SOMELIB libgcc)
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

add_executable(kurs3chm main.cpp)
