cmake_minimum_required(VERSION 3.10)

project(MyFortranProject LANGUAGES Fortran)

file(GLOB_RECURSE Fortran_SOURCES src/*.f90)

add_executable(CO2Transport ${Fortran_SOURCES})
