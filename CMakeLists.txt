cmake_minimum_required(VERSION 3.15)

set(PROJECT_NAME Maxwells_equations_FDTD)
project(${PROJECT_NAME} LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    message(STATUS "OpenMP ${OpenMP_CXX_VERSION} FOUND")
else()
    message(WARNING "OpenMP NOT found")
endif()

# set(OpenMP_omp_LIBRARY "H:/Programs/Microsoft Visual Studio/2022/Community/VC/Tools/Llvm/x64/lib/libomp.lib")
# include_directories("H:/Programs/Microsoft Visual Studio/2022/Community/VC/Tools/Llvm/x64/lib/clang/17/include")
# message(STATUS "OpenMP_omp_LIBRARY_ABOBA: ${OpenMP_omp_LIBRARY}")


# Определение путей к директории с файлами, содержащими данные для графиков от КАЖДОГО ПРОЦЕССА
set(PATH_TO_CONVERGENCE_DATA "${CMAKE_SOURCE_DIR}/input_for_graphs/convergence.csv") # Это пока хз, оставлю как есть
set(PATH_TO_CALCULATED_DATA_DIRECTORY "${CMAKE_SOURCE_DIR}/input_for_graphs/my_data/")
set(PATH_TO_ANALYTICAL_DATA_DIRECTORY "${CMAKE_SOURCE_DIR}/input_for_graphs/analytical_data/")

# # Определение путей к файлам с данными для графиков
# set(PATH_TO_CONVERGENCE_DATA "${CMAKE_SOURCE_DIR}/input_for_graphs/convergence.csv")
# set(PATH_TO_CALCULATED_DATA "${CMAKE_SOURCE_DIR}/input_for_graphs/my_data.csv")
# set(PATH_TO_ANALYTICAL_DATA "${CMAKE_SOURCE_DIR}/input_for_graphs/analytical_data.csv")

# Передача переменных в компилятор
add_definitions(-DPATH_TO_CONVERGENCE_DATA="${PATH_TO_CONVERGENCE_DATA}")
add_definitions(-DPATH_TO_CALCULATED_DATA_DIRECTORY="${PATH_TO_CALCULATED_DATA_DIRECTORY}")
add_definitions(-DPATH_TO_ANALYTICAL_DATA_DIRECTORY="${PATH_TO_ANALYTICAL_DATA_DIRECTORY}")
# add_definitions(-DPATH_TO_CALCULATED_DATA="${PATH_TO_CALCULATED_DATA}")
# add_definitions(-DPATH_TO_ANALYTICAL_DATA="${PATH_TO_ANALYTICAL_DATA}")


set(TEST_LIB GTestLib)
set(TEST_TARGET_NAME Tests)
set(ALG_TARGET_NAME Algorithms)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
# add_compile_options(-finput-charset=UTF-8 -fexec-charset=UTF-8)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/GTestLib/headers)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/algorithms/headers)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Tests/headers)


# Подключение MPI
find_package(MPI)
if (MPI_FOUND)
    message(STATUS "MPI FOUND")
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


add_subdirectory(GTestLib)
add_subdirectory(algorithms)
add_subdirectory(tests)

# Флаги компиляции OpenMP для MSVC и других компиляторов
if(MSVC AND OpenMP_CXX_FOUND)
    target_compile_options(${ALG_TARGET_NAME} PRIVATE /openmp)
    target_compile_options(${TEST_TARGET_NAME} PRIVATE /openmp)
elseif (OpenMP_CXX_FOUND)
    target_compile_options(${ALG_TARGET_NAME} PRIVATE -fopenmp)
    target_compile_options(${TEST_TARGET_NAME} PRIVATE -fopenmp)
endif()

# # Задача для запуска
# add_custom_target(run
#     COMMAND mpiexec -np 4 $<TARGET_FILE:Tests>
#     DEPENDS Tests
#     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#     COMMENT "Running Tests with MPI"
# )
