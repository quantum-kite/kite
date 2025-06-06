cmake_minimum_required(VERSION 3.9)
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++") 
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  add_definitions(-DDEBUG=1)
endif()

project(KITE)



find_package(HDF5 COMPONENTS CXX C HL)
if(${HDF5_FOUND})
  include_directories(${HDF5_INCLUDE_DIR})
else()
  MESSAGE(FATAL_ERROR  "Couldn't find hdf5")
endif()
set(HDF5_USE_STATIC_LIBRARIES ON)
MESSAGE(STATUS "Hdf5 Library:  ${HDF5_CXX_HL_LIBRARIES}")
MESSAGE(STATUS "Hdf5Hl Library:  ${HDF5_CXX_LIBRARIES}")
MESSAGE(STATUS "Compiler:  ${CMAKE_CXX_COMPILER}")
include_directories(~/include/)


find_package(Eigen3)
if(${EIGEN3_FOUND})
  include_directories(${EIGEN3_INCLUDE_DIR})
else()
  include_directories(~/include/eigen3)
  MESSAGE(STATUS "Couldn't find Eigen3. Looking in ~/include")
endif()


find_package(OpenMP)
MESSAGE(STATUS "OPENMP:  ${OPENMP_FOUND}")      
if (OPENMP_FOUND)
    set (CORRECT_CODING_FLAGS "-Wall ")
    set (CMAKE_C_FLAGS "${CORRECT_CODING_FLAGS} -g -O3 ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CORRECT_CODING_FLAGS} -g -O3 ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS_DEBUG "-g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()


set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
set(SYSTEM_NAME ${CMAKE_HOST_SYSTEM_NAME})
cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
set(MACHINE_NAME ${_host_name})
#MESSAGE(STATUS "DIR" ${PROJECT_BINARY_DIR})
include_directories(${PROJECT_BINARY_DIR}/Src)

FILE(GLOB_RECURSE SRCFILES Src/*.cpp)
#FILE(GLOB SRCFILES Src/*.cpp)
#FILE(GLOB SRCFILES2 Src/Spectral/*.cpp)
#FILE(GLOB SRCFILES3 Src/Tools/*.cpp)
#FILE(GLOB SRCFILES4 Src/CondDC/*.cpp)
#FILE(GLOB SRCFILES5 Src/OptCond_1order/*.cpp)
#FILE(GLOB SRCFILES6 Src/OptCond_2order/*.cpp)
add_executable(KITE-tools ${SRCFILES})# ${SRCFILES2} ${SRCFILES3} ${SRCFILES4} ${SRCFILES5} ${SRCFILES6})
string(TIMESTAMP TODAY "(yymmdd) %Y-%m-%d %H:%M:%S")
configure_file (
  "Src/compiletime_info.h.in"
  "Src/compiletime_info.h"
  )

include_directories(${Src})
target_link_libraries(KITE-tools ${HDF5_CXX_LIBRARIES} )
install (TARGETS KITE-tools DESTINATION bin)
