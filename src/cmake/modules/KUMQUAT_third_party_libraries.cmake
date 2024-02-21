# This files manage third party libraries required by GUDHI

find_package(Boost 1.82.0 REQUIRED)# COMPONENTS filesystem unit_test_framework program_options)

# Boost_FOUND is not reliable
if(NOT Boost_VERSION)
  message(FATAL_ERROR "NOTICE: This program requires Boost and will not be compiled.")
endif(NOT Boost_VERSION)
include_directories(${Boost_INCLUDE_DIRS})
message(STATUS "boost include dirs:" ${Boost_INCLUDE_DIRS})
message(STATUS "boost library dirs:" ${Boost_LIBRARY_DIRS})

find_package(GMP REQUIRED)
if(GMP_FOUND)
  INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR})
  message(STATUS " -- gmp include dirs:" ${GMP_INCLUDE_DIR})
  message(STATUS " -- gmp library dirs:" ${GMP_LIBRARY_DIR})
  find_package(GMPXX)
  if(GMPXX_FOUND)
    INCLUDE_DIRECTORIES(${GMPXX_INCLUDE_DIR})
    message(STATUS " -- gmpxx include dirs:" ${GMPXX_INCLUDE_DIR})
    message(STATUS " -- gmpxx library dirs:" ${GMPXX_LIBRARY_DIR})
  endif()
endif()

# from windows vcpkg eigen 3.4.0#2 : build fails with
# error C2440: '<function-style-cast>': cannot convert from 'Eigen::EigenBase<Derived>::Index' to '__gmp_expr<mpq_t,mpq_t>'
# cf. https://gitlab.com/libeigen/eigen/-/issues/2476
# Workaround is to compile with '-DEIGEN_DEFAULT_DENSE_INDEX_TYPE=int'
if (FORCE_EIGEN_DEFAULT_DENSE_INDEX_TYPE_TO_INT)
  message("++ User explicit demand to force EIGEN_DEFAULT_DENSE_INDEX_TYPE to int")
  add_definitions(-DEIGEN_DEFAULT_DENSE_INDEX_TYPE=int)
endif()

option(WITH_GUDHI_USE_TBB "Build with Intel TBB parallelization" ON)

# Find TBB package for parallel sort - not mandatory, just optional.
#if(WITH_KUMQUAT_USE_TBB)
find_package(TBB CONFIG)
if(TARGET TBB::tbb)
  # Specific windows case with its debug/release management
  if(CMAKE_BUILD_TYPE MATCHES Debug)
    get_target_property(TBB_LIBRARY TBB::tbb IMPORTED_LOCATION_DEBUG)
    get_target_property(TBB_MALLOC_LIBRARY TBB::tbbmalloc IMPORTED_LOCATION_DEBUG)
  else()
    get_target_property(TBB_LIBRARY TBB::tbb IMPORTED_LOCATION_RELEASE)
    get_target_property(TBB_MALLOC_LIBRARY TBB::tbbmalloc IMPORTED_LOCATION_RELEASE)
  endif()
  # Generic case
  if (NOT TBB_LIBRARY)
    get_target_property(TBB_LIBRARY TBB::tbb LOCATION)
    get_target_property(TBB_MALLOC_LIBRARY TBB::tbbmalloc LOCATION)
  endif()
  # TBB Error management
  if (TBB_VERSION VERSION_LESS 2019.0.11007)
    # TBBTargets.cmake was introduced in 2019.7, so this case should not happen
    # cf. https://github.com/oneapi-src/oneTBB/blob/2019_U7/CHANGES
    message(WARNING "++ TBB found but version ${TBB_VERSION} is too old - KUMQUAT cannot compile with TBB")
  else()
    if (NOT TBB_LIBRARY)
      message(WARNING "++ TBB found but not ${TBB_LIBRARY} - KUMQUAT cannot compile with TBB")
    else()
      # A correct version of TBB was found
      get_target_property(TBB_INCLUDE_DIRS TBB::tbb INTERFACE_INCLUDE_DIRECTORIES)
      get_filename_component(TBB_LIBRARY_DIRS ${TBB_LIBRARY} DIRECTORY)
      message("++ TBB version ${TBB_VERSION} found in ${TBB_LIBRARY_DIRS} - includes in ${TBB_INCLUDE_DIRS}")
      add_definitions(-DGUDHI_USE_TBB)
      if(MSVC)
        # cf. https://github.com/oneapi-src/oneTBB/issues/573
        add_definitions(-DNOMINMAX)
      endif()
    endif()
  endif()
endif()
#endif()

find_package(Eigen3 3.1.0)
if (EIGEN3_FOUND)
  include( ${EIGEN3_USE_FILE} )
  set(CGAL_WITH_EIGEN3_VERSION ${CGAL_VERSION})
endif (EIGEN3_FOUND)

# Required programs for unitary tests purpose
FIND_PROGRAM( GCOVR_PATH gcovr )
if (GCOVR_PATH)
  message("gcovr found in ${GCOVR_PATH}")
endif()
FIND_PROGRAM( GPROF_PATH gprof )
if (GPROF_PATH)
  message("gprof found in ${GPROF_PATH}")
endif()
FIND_PROGRAM( DIFF_PATH diff )
if (DIFF_PATH)
  message("diff found in ${DIFF_PATH}")
endif()
FIND_PROGRAM( GNUPLOT_PATH gnuplot )
if (GNUPLOT_PATH)
  message("gnuplot found in ${GNUPLOT_PATH}")
endif()

# BOOST ISSUE result_of vs C++11
add_definitions(-DBOOST_RESULT_OF_USE_DECLTYPE)
# BOOST ISSUE with Libraries name resolution under Windows
add_definitions(-DBOOST_ALL_NO_LIB)
# problem with Visual Studio link on Boost program_options
add_definitions( -DBOOST_ALL_DYN_LINK )
# problem on Mac with boost_system and boost_thread
add_definitions( -DBOOST_SYSTEM_NO_DEPRECATED )