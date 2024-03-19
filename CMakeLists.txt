cmake_minimum_required(VERSION 3.0.0)

set(PROJECT_NAME Maxwells_equations_FDTD)
project(${PROJECT_NAME} LANGUAGES CXX C)

find_package(OpenMP)


set(testLib GTestLib)
set(test Tests)
set(alg Algorithms)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/GTestLib/headers)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/algorithms/headers)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Tests/headers)

set(MPI_INCLUDE_DIR H:/Microsoft_MPI/SDK/Include)
set(MPI_LIB_DIR 	H:/Microsoft_MPI/SDK/Lib/x64)
set(MPI_LIBS 		H:/Microsoft_MPI/SDK/Lib/x64/msmpi.lib)

include_directories(${MPI_INCLUDE_DIR})


add_subdirectory(GTestLib)
add_subdirectory(algorithms)
add_subdirectory(tests)

#if(OpenMP_CXX_FOUND)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /openmp")
#    target_link_libraries(${test} OpenMP)
#endif()
