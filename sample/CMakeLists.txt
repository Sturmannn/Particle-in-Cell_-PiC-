file(GLOB hdrs "*.h*")
file(GLOB srcs "*.cpp")

add_executable(${FDTD_SAMPLE} ${srcs} ${hdrs})

target_include_directories(${FDTD_SAMPLE} PUBLIC 
    ${CMAKE_SOURCE_DIR}/include/core
    ${CMAKE_SOURCE_DIR}/include/parallel
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(${FDTD_SAMPLE} PUBLIC  ${FDTD_LIB} MPI::MPI_CXX  OpenMP::OpenMP_CXX)
else()
    target_link_libraries(${FDTD_SAMPLE} PUBLIC  ${FDTD_LIB} MPI::MPI_CXX)
endif()

