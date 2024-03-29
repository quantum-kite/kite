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
  set(HDF5_USE_STATIC_LIBRARIES ON)
  MESSAGE(STATUS "Hdf5 Library:  ${HDF5_CXX_HL_LIBRARIES}")
  MESSAGE(STATUS "Hdf5Hl Library:  ${HDF5_CXX_LIBRARIES}")
else()
  MESSAGE(STATUS "Couldn't find HDF5")
endif()
MESSAGE(STATUS "Compiler:  ${CMAKE_CXX_COMPILER}")
include_directories(~/include/)

find_package(Eigen3)
if(${EIGEN3_FOUND})
  include_directories(${EIGEN3_INCLUDE_DIR})
else()
  include_directories(~/include/eigen3)
  MESSAGE(STATUS "Couldn't find Eigen3. Looking in ~/include")
endif()

set (CORRECT_CODING_FLAGS "")
find_package(OpenMP)
MESSAGE(STATUS "OPENMP:  ${OPENMP_FOUND}")      
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CORRECT_CODING_FLAGS} ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CORRECT_CODING_FLAGS} -O3 ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CORRECT_CODING_FLAGS} ${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CORRECT_CODING_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "${CORRECT_CODING_FLAGS} -g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()






FILE(GLOB_RECURSE SRCFILES Src/*.cpp)
MESSAGE(STATUS "Files to compile" ${SRCFILES})
add_executable(KITEx ${SRCFILES})
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

set(compile_wp 1)
if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0.0)
    MESSAGE( STATUS "Cannot compile gaussian wavepacket")
    set(compile_wp 0)
endif()
add_definitions(-DCOMPILE_WAVEPACKET=${compile_wp})


target_link_libraries(KITEx ${HDF5_CXX_LIBRARIES})
target_include_directories(KITEx PRIVATE "Src/"
"Src/Tools/"
"Src/Lattice/"
"Src/Hamiltonian/"
"Src/Vector/"
"Src/Simulation/"
)

install (TARGETS KITEx DESTINATION bin)



set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
set(SYSTEM_NAME ${CMAKE_HOST_SYSTEM_NAME})
cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
set(MACHINE_NAME ${_host_name})

include_directories(${PROJECT_BINARY_DIR}/tools/Src)
FILE(GLOB_RECURSE SRCFILES tools/Src/*.cpp)
add_executable(KITE-tools ${SRCFILES})
string(TIMESTAMP TODAY "(yymmdd) %Y-%m-%d %H:%M:%S")
configure_file (
  "tools/Src/compiletime_info.h.in"
  "tools/Src/compiletime_info.h"
  )

include_directories(${Src})
target_link_libraries(KITE-tools ${HDF5_CXX_LIBRARIES} )
install (TARGETS KITE-tools DESTINATION bin)
