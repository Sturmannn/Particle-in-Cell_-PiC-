file(GLOB core_hdrs "core/*.h*")
file(GLOB parallel_hdrs "parallel/*.h*")
file(GLOB shared_hdrs "shared/*.h*")

file(GLOB core_srcs "core/*.cpp")
file(GLOB parallel_srcs "parallel/*.cpp")
file(GLOB shared_srcs "shared/*.cpp")

add_library(${FDTD_LIB} STATIC ${core_hdrs} ${parallel_hdrs} ${shared_hdrs} ${core_srcs} ${parallel_srcs} ${shared_srcs})

target_include_directories(${FDTD_LIB} PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/core
    ${CMAKE_CURRENT_SOURCE_DIR}/parallel
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(${FDTD_LIB} PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX)
    # target_compile_options(${FDTD_LIB} PRIVATE ${CMAKE_CXX_FLAGS} -fopenmp)
else()
    target_link_libraries(${FDTD_LIB} PUBLIC MPI::MPI_CXX)
    target_compile_options(${FDTD_LIB} PRIVATE ${CMAKE_CXX_FLAGS})
endif()
