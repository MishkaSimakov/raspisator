add_library(raspisator INTERFACE)

target_include_directories(raspisator INTERFACE .)
target_link_libraries(raspisator INTERFACE fmt::fmt)