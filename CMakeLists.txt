cmake_minimum_required(VERSION 3.15)

set(PROJECT_NAME Maxwells_equations_FDTD)
project(${PROJECT_NAME} LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -qopenmp -g")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3 -fopt-info-vec -fopt-info-vec-missed")


# Определение путей к директории с файлами, содержащими данные для графиков от КАЖДОГО ПРОЦЕССА
set(PATH_TO_CONVERGENCE_DATA "${CMAKE_SOURCE_DIR}/input_for_graphs/convergence.csv") # Это пока хз, оставлю как есть
set(PATH_TO_CALCULATED_DATA_DIRECTORY "${CMAKE_SOURCE_DIR}/input_for_graphs/my_data/")
set(PATH_TO_ANALYTICAL_DATA_DIRECTORY "${CMAKE_SOURCE_DIR}/input_for_graphs/analytical_data/")
set(PATH_TO_MEASUREMENTS_FILE "${CMAKE_SOURCE_DIR}/measurements.csv")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

# Передача переменных в компилятор
add_definitions(-DPATH_TO_CONVERGENCE_DATA="${PATH_TO_CONVERGENCE_DATA}")
add_definitions(-DPATH_TO_CALCULATED_DATA_DIRECTORY="${PATH_TO_CALCULATED_DATA_DIRECTORY}")
add_definitions(-DPATH_TO_ANALYTICAL_DATA_DIRECTORY="${PATH_TO_ANALYTICAL_DATA_DIRECTORY}")
add_definitions(-DPATH_TO_MEASUREMENTS_FILE="${PATH_TO_MEASUREMENTS_FILE}")


# set(TEST_LIB GTestLib)
set(TEST_TARGET_NAME Tests)
set(ALG_TARGET_NAME Algorithms)

# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/googletest/include)
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/GTestLib/headers)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/algorithms/headers)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Tests/headers)

# find_package(OpenMP REQUIRED)
# if(OpenMP_CXX_FOUND)
#     message(STATUS "OpenMP C++ FOUND - ${OpenMP_CXX_VERSION}")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# else()
#     message(WARNING "OpenMP NOT found")
# endif()

# Подключение MPI
find_package(MPI)
if (MPI_CXX_FOUND)
    message(STATUS "MPI C++ FOUND - ${MPI_CXX_VERSION}")
    # set(MPI_INCLUDE_DIR ${MPI_INCLUDE_PATH})
else()
    message(WARNING "MPI NOT found")

    set(MPI_INCLUDE_DIR H:/Microsoft_MPI/SDK/Include)
    set(MPI_LIB_DIR 	H:/Microsoft_MPI/SDK/Lib/x64)
    set(MPI_LIBS 		H:/Microsoft_MPI/SDK/Lib/x64/msmpi.lib)
    
    include_directories(${MPI_INCLUDE_DIR})
    
    message(STATUS "MPI_INCLUDE_DIR: ${MPI_INCLUDE_DIR}")
    message(STATUS "MPI_LIB_DIR: ${MPI_LIB_DIR}")
    message(STATUS "MPI_LIBS: ${MPI_LIBS}")
endif()


# add_subdirectory(GTestLib)

option(BUILD_GMOCK "Build gmock" OFF)
add_subdirectory(3rdparty/googletest)

add_subdirectory(algorithms)
add_subdirectory(tests)

# Флаги компиляции OpenMP для MSVC и других компиляторов
# if(MSVC AND OpenMP_CXX_FOUND)
#     target_compile_options(${ALG_TARGET_NAME} PRIVATE /openmp)
#     target_compile_options(${TEST_TARGET_NAME} PRIVATE /openmp)
# elseif (OpenMP_CXX_FOUND)
#     target_compile_options(${ALG_TARGET_NAME} PRIVATE -fopenmp)
#     target_compile_options(${TEST_TARGET_NAME} PRIVATE -fopenmp)
# endif()

# Проверка векторизации
if (CMAKE_CXX_COMPILER_ID MATCHES  "Clang")
    message(STATUS "Clang compiler")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -g")
    # add_compile_options(-Rpass=loop-vectorize)
    # add_compile_options(-Rpass-missed=loop-vectorize)
    # add_compile_options(-Rpass-analysis=loop-vectorize)
elseif (CMAKE_CXX_COMPILER_ID MATCHES  "GNU")
# add_compile_options(-fopt-info-vec-optimized)
    message(STATUS "GNU compiler")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -g -fopt-info-vec -fopt-info-vec-missed")
elseif (CMAKE_CXX_COMPILER_ID MATCHES  "Intel")
    message(STATUS "Intel compiler")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3 -g -qopt-report=5")
elseif (CMAKE_CXX_COMPILER_ID MATCHES  "MSVC")
    message(STATUS "MSVC compiler")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O3 /openmp")
endif()

message(STATUS "Using compiler: " ${CMAKE_CXX_COMPILER_ID})
message(STATUS "CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})
