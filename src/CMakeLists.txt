add_library(raspisator linear/Matrix.cpp)

target_include_directories(raspisator PUBLIC .)
target_link_libraries(raspisator PUBLIC fmt::fmt)