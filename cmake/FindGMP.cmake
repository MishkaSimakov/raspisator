set(GMP_PREFIX "" CACHE PATH "GMP path prefix")

find_path(GMP_INCLUDE_DIR gmp.h gmpxx.h
        PATHS ${GMP_PREFIX}/include /usr/include /usr/local/include)

find_library(GMP_LIBRARY NAMES gmp libgmp
        PATHS ${GMP_PREFIX}/lib /usr/lib /usr/local/lib)

find_library(GMPXX_LIBRARY NAMES gmpxx libgmpxx
        PATHS ${GMP_PREFIX}/lib /usr/lib /usr/local/lib)

if (GMP_INCLUDE_DIR AND GMP_LIBRARY AND GMPXX_LIBRARY)
    set(GMP_FOUND TRUE)
endif ()

if (GMP_FOUND)
    if (NOT GMP_FIND_QUIETLY)
        message(STATUS "Found GMP: ${GMP_LIBRARY} ${GMPXX_LIBRARY}")
    endif ()
elseif (GMP_FOUND)
    if (GMP_FIND_REQUIRED)
        message(FATAL_ERROR "Could not find GMP")
    endif ()
endif ()