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
