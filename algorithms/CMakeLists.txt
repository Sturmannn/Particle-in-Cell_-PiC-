file(GLOB hdrs "headers/*.h*")
file(GLOB srcs "sources/*.cpp")

add_library(${ALG_TARGET_NAME} STATIC ${srcs} ${hdrs})

target_include_directories(${ALG_TARGET_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/headers)

if(OpenMP_CXX_FOUND)
    target_link_libraries(${ALG_TARGET_NAME} PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX)
else()
    target_link_libraries(${ALG_TARGET_NAME} PUBLIC MPI::MPI_CXX)
endif()

