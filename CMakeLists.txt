cmake_minimum_required(VERSION 3.15)

set(PROJECT_NAME Maxwells_equations_FDTD)
project(${PROJECT_NAME} LANGUAGES CXX C)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
message(STATUS "OpenMP ${OpenMP_CXX_VERSION} FOUND")
else()
message(WARNING "OpenMP NOT found")
endif()

# Определение путей к файлам с данными для графиков
set(PATH_TO_CONVERGENCE_DATA "${CMAKE_SOURCE_DIR}/input_for_graphs/convergence.csv")
set(PATH_TO_CALCULATED_DATA "${CMAKE_SOURCE_DIR}/input_for_graphs/my_data.csv")
set(PATH_TO_ANALYTICAL_DATA "${CMAKE_SOURCE_DIR}/input_for_graphs/analytical_data.csv")

# Передача переменных в компилятор
add_definitions(-DPATH_TO_CONVERGENCE_DATA="${PATH_TO_CONVERGENCE_DATA}")
add_definitions(-DPATH_TO_CALCULATED_DATA="${PATH_TO_CALCULATED_DATA}")
add_definitions(-DPATH_TO_ANALYTICAL_DATA="${PATH_TO_ANALYTICAL_DATA}")

set(testLib GTestLib)
set(test Tests)
set(alg Algorithms)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
# add_compile_options(-finput-charset=UTF-8 -fexec-charset=UTF-8)

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
