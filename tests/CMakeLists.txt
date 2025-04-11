file(GLOB hdrs "headers/*.h*") #Добавляем в переменную hdrs все файлы с расширением .h
file(GLOB srcs "sources/*.cpp") #Добавляем в переменную srcs все файлы с расширением .cpp

add_executable(${TEST_TARGET_NAME} ${srcs} ${hdrs})

# Добавляем путь к заголовочным файлам
target_include_directories(${TEST_TARGET_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/headers)
target_include_directories(${TEST_TARGET_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/3rdparty/googletest/googletest/include)


#Подключаем к проекту с именем ${TEST_TARGET_NAME} библиотеку с именем ${testLib}
if(OpenMP_CXX_FOUND)
    # target_link_libraries(${TEST_TARGET_NAME} PUBLIC ${TEST_LIB} ${ALG_TARGET_NAME} ${MPI_LIBS}  OpenMP::OpenMP_CXX)
    target_link_libraries(${TEST_TARGET_NAME} PUBLIC  ${ALG_TARGET_NAME} MPI::MPI_CXX  OpenMP::OpenMP_CXX GTest::gtest)
else()
    # target_link_libraries(${TEST_TARGET_NAME} PUBLIC ${TEST_LIB} ${ALG_TARGET_NAME} ${MPI_LIBS})
    target_link_libraries(${TEST_TARGET_NAME} PUBLIC  ${ALG_TARGET_NAME} MPI::MPI_CXX GTest::gtest)
endif()


# target_link_libraries(${TEST_TARGET_NAME} PUBLIC ${TEST_LIB} ${ALG_TARGET_NAME} ${MPI_LIBS} ${OpenMP_omp_LIBRARY})
