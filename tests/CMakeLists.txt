file(GLOB hdrs "headers/*.h*")
file(GLOB srcs "sources/*.cpp")

add_executable(${TEST_TARGET_NAME} ${srcs} ${hdrs})

target_include_directories(${TEST_TARGET_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/headers)
target_include_directories(${TEST_TARGET_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/3rdparty/googletest/googletest/include)

if(OpenMP_CXX_FOUND)
    target_link_libraries(${TEST_TARGET_NAME} PUBLIC  ${ALG_TARGET_NAME} MPI::MPI_CXX  OpenMP::OpenMP_CXX GTest::gtest)
else()
    target_link_libraries(${TEST_TARGET_NAME} PUBLIC  ${ALG_TARGET_NAME} MPI::MPI_CXX GTest::gtest)
endif()

