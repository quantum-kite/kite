cmake_minimum_required (VERSION 3.11.4)
project (Post-processing)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_STANDARD 17)


find_package(HDF5 COMPONENTS CXX C HL)
set(HDF5_USE_STATIC_LIBRARIES ON)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(OpenMP)
    MESSAGE(STATUS "OPENMP:  ${OPENMP_FOUND}")      
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")
endif()


#add_compile_options(-ffast-math -O2 -DEIGEN_DONT_PARALLELIZE)
add_compile_options(-Wall -ffast-math -march=native -O2 -fdiagnostics-color=always -DEIGEN_DONT_PARALLELIZE)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

MESSAGE(STATUS "Hdf5HL Library:  ${HDF5_CXX_HL_LIBRARIES}")
MESSAGE(STATUS "Hdf5 Library:  ${HDF5_CXX_LIBRARIES}")
MESSAGE(STATUS "Compiler:  ${CMAKE_CXX_COMPILER}")

FILE(GLOB sources "src/*.cpp")
FILE(GLOB sources2 "src/cond_2order/*.cpp")
FILE(GLOB includes "src/")
FILE(GLOB includes2 "src/cond_2order/")
#FILE(GLOB includes "src/*.hpp")
#FILE(GLOB includes2 "src/cond_2order/*.hpp")

add_executable(KITE-tools ${sources} ${sources2})
include_directories(${includes} ${includes2})
#add_executable(KITE-tools ${includes} ${includes2} ${sources} ${sources2})

target_link_libraries(KITE-tools ${HDF5_CXX_LIBRARIES})
