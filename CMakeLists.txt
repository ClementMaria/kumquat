cmake_minimum_required(VERSION 3.8)

project(KUMQUAT)

set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/src/cmake/modules/")

find_package(Boost 1.71.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
message(STATUS " -- boost include dirs:" ${Boost_INCLUDE_DIRS})
message(STATUS " -- boost library dirs:" ${Boost_LIBRARY_DIRS})

find_package(GMP REQUIRED)
if(GMP_FOUND)
  include_directories(${GMP_INCLUDE_DIR})
  message(STATUS " -- gmp include dirs:" ${GMP_INCLUDE_DIR})
  message(STATUS " -- gmp library dirs:" ${GMP_LIBRARY_DIR})
  find_package(GMPXX REQUIRED)
  if(GMPXX_FOUND)
    include_directories(${GMPXX_INCLUDE_DIR})
    message(STATUS " -- gmpxx include dirs:" ${GMPXX_INCLUDE_DIR})
    message(STATUS " -- gmpxx library dirs:" ${GMPXX_LIBRARY_DIR})
  endif()
endif()

find_package(Eigen3 REQUIRED)
include_directories(${Eigen_INCLUDE_DIRS})
# link_directories(${Eigen_INCLUDE_DIRS})

include(KUMQUAT_modules)

add_kumquat_module(algebra)
add_kumquat_module(matrix)

# Include module CMake subdirectories
# KUMQUAT_SUB_DIRECTORIES is managed in CMAKE_MODULE_PATH/KUMQUAT_modules.cmake
foreach(KUMQUAT_MODULE ${KUMQUAT_MODULES})
  foreach(KUMQUAT_SUB_DIRECTORY ${KUMQUAT_SUB_DIRECTORIES})
    if(EXISTS ${CMAKE_SOURCE_DIR}/src/${KUMQUAT_MODULE}/${KUMQUAT_SUB_DIRECTORY}/CMakeLists.txt)
      add_subdirectory(src/${KUMQUAT_MODULE}/${KUMQUAT_SUB_DIRECTORY}/)
    endif()
  endforeach()
endforeach()