cmake_minimum_required(VERSION 3.9)

project(KUMQUAT)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/src/cmake/modules/")

set(CMAKE_CXX_STANDARD 23)

# Reset cache
set(KUMQUAT_MODULES "" CACHE INTERNAL "KUMQUAT_MODULES")
set(KUMQUAT_MISSING_MODULES "" CACHE INTERNAL "KUMQUAT_MISSING_MODULES")

include(KUMQUAT_third_party_libraries NO_POLICY_SCOPE)

# add_compile_options(-Wall -Wextra -Wpedantic -pthread -lmpc -lmpfr -lgmp)

include(KUMQUAT_modules)

add_kumquat_module(algebraic_structures)
add_kumquat_module(algebraic_elements)
add_kumquat_module(matrix)
add_kumquat_module(utils)
add_kumquat_module(random_process)
add_kumquat_module(braid)
add_kumquat_module(quantum_topology)

# Include module CMake subdirectories
# KUMQUAT_SUB_DIRECTORIES is managed in CMAKE_MODULE_PATH/KUMQUAT_modules.cmake
foreach(KUMQUAT_MODULE ${KUMQUAT_MODULES})
  foreach(KUMQUAT_SUB_DIRECTORY ${KUMQUAT_SUB_DIRECTORIES})
    if(EXISTS ${CMAKE_SOURCE_DIR}/src/${KUMQUAT_MODULE}/${KUMQUAT_SUB_DIRECTORY}/CMakeLists.txt)
      add_subdirectory(src/${KUMQUAT_MODULE}/${KUMQUAT_SUB_DIRECTORY}/)
    endif()
  endforeach()
endforeach()

#include(KUMQUAT_doxygen_target)

message("++ KUMQUAT_MODULES list is:\"${KUMQUAT_MODULES}\"")

message("++-- ${Boost_LIBRARIES}   |   ${GMPXX_LIBRARIES}   |   ${GMP_LIBRARIES}")