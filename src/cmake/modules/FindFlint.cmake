# include(LibFindMacros)

function (libfind_library libname pkg)
    string(TOUPPER ${pkg} PKG)
    string(TOUPPER ${libname} LIBNAME)

    find_library(${LIBNAME}_LIBRARY
        NAMES
            ${libname} ${ARGN}
    )

    if (NOT TARGET ${libname})
        add_library(${libname} UNKNOWN IMPORTED)
        set_property(TARGET ${libname} PROPERTY IMPORTED_LOCATION ${${LIBNAME}_LIBRARY})
    endif()
endfunction()

function (libfind_include HEADER pkg)
    string(TOUPPER ${pkg} PKG)

    find_path(${PKG}_INCLUDE_DIR
        NAMES
            ${HEADER} ${ARGN}
    )
endfunction()

libfind_include(flint/flint.h flint)
libfind_library(flint flint)

set(FLINT_LIBRARIES ${FLINT_LIBRARY})
# Workaround for https://github.com/fredrik-johansson/arb/issues/24
set(FLINT_INCLUDE_DIRS ${FLINT_INCLUDE_DIR} ${FLINT_INCLUDE_DIR}/flint)

set(FLINT_TARGETS flint)

file(READ "${FLINT_INCLUDE_DIR}/flint/flint.h" FLINT_H)
string(REGEX MATCH "#define FLINT_VERSION \"([0-9.]*)" _ ${FLINT_H})
set(FLINT_VERSION ${CMAKE_MATCH_1})
string(REPLACE "." ";" FLINT_VERSION_LIST ${FLINT_VERSION})
list(GET FLINT_VERSION_LIST 0 FLINT_VERSION_MAJOR)
list(GET FLINT_VERSION_LIST 1 FLINT_VERSION_MINOR)
list(GET FLINT_VERSION_LIST 2 FLINT_VERSION_PATCH)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FLINT DEFAULT_MSG FLINT_LIBRARIES
  FLINT_INCLUDE_DIRS FLINT_VERSION)

mark_as_advanced(FLINT_INCLUDE_DIR FLINT_LIBRARY FLINT_H FLINT_VERSION_LIST)