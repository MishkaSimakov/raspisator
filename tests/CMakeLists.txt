enable_testing()

add_executable(tests)

target_sources(tests PUBLIC BasicTests.cpp)
target_link_libraries(tests GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(tests)