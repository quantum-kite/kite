cmake_minimum_required(VERSION 2.8)
project(KITEx)
set(CMAKE_C_COMPILER "/opt/local/bin/gcc-mp-8")
SET(CMAKE_CXX_COMPILER "/opt/local/bin/g++-mp-8") 
set(CMAKE_CXX_STANDARD 17)

FILE(GLOB SRCFILES *.cpp)
add_executable(KITEx ${SRCFILES})

find_package(HDF5 COMPONENTS CXX C HL)
set(HDF5_USE_STATIC_LIBRARIES ON)
MESSAGE(STATUS "Hdf5 Library:  ${HDF5_CXX_HL_LIBRARIES}")
MESSAGE(STATUS "Hdf5Hl Library:  ${HDF5_CXX_LIBRARIES}")
MESSAGE(STATUS "Compiler:  ${CMAKE_CXX_COMPILER}")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(OpenMP)
    MESSAGE(STATUS "OPENMP:  ${OPENMP_FOUND}")      
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")
endif()


add_compile_options(-Wall -ffast-math -O2 -fdiagnostics-color=always -DEIGEN_DONT_PARALLELIZE)


include_directories(${Src})
target_link_libraries(KITEx ${HDF5_CXX_LIBRARIES} )




