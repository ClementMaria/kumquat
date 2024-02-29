# This files manage compilation flags required by GUDHI

include(TestCXXAcceptsFlag)

# add a compiler flag only if it is accepted
macro(add_cxx_compiler_flag _flag)
  string(REPLACE "-" "_" "/" _flag_var ${_flag})
  check_cxx_accepts_flag("${_flag}" CXX_COMPILER_${_flag_var}_OK)
  if(CXX_COMPILER_${_flag_var}_OK)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_flag}")
  endif()
endmacro()

set (CMAKE_CXX_STANDARD 20)

try_compile(CXX20_VECTOR_TEMPLATE_DEDUCTION
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_LIST_DIR}/vector_template_deduction.cpp
  CMAKE_FLAGS -DCMAKE_CXX_STANDARD=20
  OUTPUT_VARIABLE OUTPUT
  )

if (NOT CXX20_VECTOR_TEMPLATE_DEDUCTION)
  message(FATAL_ERROR "Your compiler does not support c++20, this is mandatory to compile the KUMQUAT library. ${CXX20_VECTOR_TEMPLATE_DEDUCTION} - ${OUTPUT}")
endif()

enable_testing()

if(MSVC)
  add_cxx_compiler_flag("/W3")
else()
  add_cxx_compiler_flag("-Wall")
  # Only for dev version
  if(PROJECT_NAME STREQUAL "KUMQUATdev")
    add_cxx_compiler_flag("-pedantic")
  endif()
endif()

if (DEBUG_TRACES)
  # For programs to be more verbose
  message(STATUS "DEBUG_TRACES are activated")
  add_definitions(-DDEBUG_TRACES)
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
  message("++ Debug compilation flags are: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
else()
  message("++ Release compilation flags are: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

message("++ Compilation flags are: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
