cmake_minimum_required(VERSION 3.10)

project(MyFortranProject LANGUAGES Fortran)

file(GLOB_RECURSE Fortran_SOURCES src/*.f90)

# Find LAPACK and BLAS libraries
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

add_executable(CO2phaseReactions ${Fortran_SOURCES})

# Link against LAPACK and BLAS
target_link_libraries(CO2phaseReactions ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
