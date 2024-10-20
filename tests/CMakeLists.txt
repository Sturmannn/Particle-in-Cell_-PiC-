file(GLOB hdrs "headers/*.h*") #Добавляем в переменную hdrs все файлы с расширением .h
file(GLOB srcs "sources/*.cpp") #Добавляем в переменную srcs все файлы с расширением .cpp

add_executable(${TEST_TARGET_NAME} ${srcs} ${hdrs})

#Подключаем к проекту с именем ${TEST_TARGET_NAME} библиотеку с именем ${testLib}
if(OpenMP_CXX_FOUND)
    target_link_libraries(${TEST_TARGET_NAME} PUBLIC ${TEST_LIB} ${ALG_TARGET_NAME} ${MPI_LIBS}  OpenMP::OpenMP_CXX)
else()
    target_link_libraries(${TEST_TARGET_NAME} PUBLIC ${TEST_LIB} ${ALG_TARGET_NAME} ${MPI_LIBS})
endif()


# target_link_libraries(${TEST_TARGET_NAME} PUBLIC ${TEST_LIB} ${ALG_TARGET_NAME} ${MPI_LIBS} ${OpenMP_omp_LIBRARY})
