file(GLOB hdrs "headers/*.h*") #Добавляем в переменную hdrs все файлы с расширением .h
file(GLOB srcs "sources/*.cpp") #Добавляем в переменную srcs все файлы с расширением .cpp

add_library(${ALG_TARGET_NAME} STATIC ${srcs} ${hdrs})

if(OpenMP_CXX_FOUND)
    # target_link_libraries(${ALG_TARGET_NAME} PUBLIC ${MPI_LIBS} OpenMP::OpenMP_CXX)
    target_link_libraries(${ALG_TARGET_NAME} PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX)
else()
    # target_link_libraries(${ALG_TARGET_NAME} PUBLIC ${MPI_LIBS})
    target_link_libraries(${ALG_TARGET_NAME} PUBLIC MPI::MPI_CXX)
endif()

# target_link_libraries(${ALG_TARGET_NAME} PUBLIC ${OpenMP_omp_LIBRARY})
