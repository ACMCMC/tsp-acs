cmake_minimum_required(VERSION 2.6.0)

# here we specify that the project is C language only, so the default
# C compiler on the system will be used
project(tsp)

file(GLOB SOURCES
    headers/*.h
    implementations/*.cpp
)

add_executable(tsp ${SOURCES})

enable_testing()
#add_subdirectory(tests)