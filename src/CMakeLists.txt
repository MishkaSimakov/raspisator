add_library(raspisator utils/Drawing.cpp linear/Matrix.cpp linear/SimplexMethod.cpp)

target_include_directories(raspisator PUBLIC .)
target_link_libraries(raspisator PUBLIC fmt::fmt)