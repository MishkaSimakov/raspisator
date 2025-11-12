include(FetchContent)

# Google Test
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/52eb8108c5bdec04579160ae17225d66034bd723.zip
)
FetchContent_MakeAvailable(googletest)


# fmt
FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt
        GIT_TAG 407c905e45ad75fc29bf0f9bb7c5c2fd3475976f
)
FetchContent_MakeAvailable(fmt)
