cmake_minimum_required(VERSION 3.0.0)

set(PROJECT_NAME Maxwells_equations_FDTD)
project(${PROJECT_NAME} LANGUAGES CXX C)

set(testLib GTestLib)
set(test Tests)
set(alg Algorithms)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/GTestLib/headers)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/algorithms/headers)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Tests/headers)


add_subdirectory(GTestLib)
add_subdirectory(algorithms)
add_subdirectory(tests)
