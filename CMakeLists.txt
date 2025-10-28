cmake_minimum_required(VERSION 3.10)

project(MyFortranProject LANGUAGES Fortran)
# Find LAPACK and BLAS libraries
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

file(GLOB_RECURSE Fortran_SOURCES src/*.f90)

set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch -U")

add_executable(CO2Transport ${Fortran_SOURCES})
target_link_libraries(CO2Transport ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
