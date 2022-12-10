cmake_minimum_required(VERSION 2.6.0)

file(GLOB TESTS
    ./*.cpp
)

add_executable(tsp_tests ${TESTS})
add_test(tsp_test tsp_tests)