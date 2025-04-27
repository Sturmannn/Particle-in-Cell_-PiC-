cmake_minimum_required(VERSION 3.15)

set(PROJECT_NAME Maxwells_equations_FDTD)
project(${PROJECT_NAME} LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

# Targets
set(TEST_TARGET_NAME Tests)
set(ALG_TARGET_NAME Algorithms)

# Определение путей к директории с файлами, содержащими данные для графиков от КАЖДОГО ПРОЦЕССА
set(PATH_TO_CONVERGENCE_DATA "${CMAKE_SOURCE_DIR}/input_for_graphs/convergence.csv") # Это пока хз, оставлю как есть
set(PATH_TO_CALCULATED_DATA_DIRECTORY "${CMAKE_SOURCE_DIR}/input_for_graphs/my_data/")
set(PATH_TO_ANALYTICAL_DATA_DIRECTORY "${CMAKE_SOURCE_DIR}/input_for_graphs/analytical_data/")
set(PATH_TO_MEASUREMENTS_FILE "${CMAKE_SOURCE_DIR}/measurements.csv")


# Передача переменных в компилятор
add_definitions(-DPATH_TO_CONVERGENCE_DATA="${PATH_TO_CONVERGENCE_DATA}")
add_definitions(-DPATH_TO_CALCULATED_DATA_DIRECTORY="${PATH_TO_CALCULATED_DATA_DIRECTORY}")
add_definitions(-DPATH_TO_ANALYTICAL_DATA_DIRECTORY="${PATH_TO_ANALYTICAL_DATA_DIRECTORY}")
add_definitions(-DPATH_TO_MEASUREMENTS_FILE="${PATH_TO_MEASUREMENTS_FILE}")

# OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP C++ FOUND - ${OpenMP_CXX_VERSION}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(WARNING "OpenMP NOT found")
endif()

# MPI
find_package(MPI REQUIRED)
if (MPI_CXX_FOUND)
    message(STATUS "MPI C++ FOUND - ${MPI_CXX_VERSION}")
else()
    message(WARNING "MPI NOT found")
endif()

# # Vectorization check
# if (CMAKE_CXX_COMPILER_ID MATCHES  "Clang")
#     message(STATUS "Clang compiler")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -g")
#     # add_compile_options(-Rpass=loop-vectorize)
#     # add_compile_options(-Rpass-missed=loop-vectorize)
#     # add_compile_options(-Rpass-analysis=loop-vectorize)
# elseif (CMAKE_CXX_COMPILER_ID MATCHES  "GNU")
# # add_compile_options(-fopt-info-vec-optimized)
#     message(STATUS "GNU compiler")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -g -fopt-info-vec -fopt-info-vec-missed")
# elseif (CMAKE_CXX_COMPILER_ID MATCHES  "Intel")
#     message(STATUS "Intel compiler")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3 -g -qopt-report=5")
# elseif (CMAKE_CXX_COMPILER_ID MATCHES  "MSVC")
#     message(STATUS "MSVC compiler")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O3 /openmp")
# endif()

option(BUILD_GMOCK "Build gmock" OFF)
add_subdirectory(3rdparty/googletest)

add_subdirectory(algorithms)
add_subdirectory(tests)


message(STATUS "Using compiler: " ${CMAKE_CXX_COMPILER_ID})
message(STATUS "CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

# ---------------------------------------------------------------------
set (FDTD_LIB "LibFDTD")
set (FDTD_SAMPLE "SampleFDTD")
add_subdirectory(include)
add_subdirectory(sample)