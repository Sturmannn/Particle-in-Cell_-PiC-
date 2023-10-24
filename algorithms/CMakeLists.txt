file(GLOB hdrs "headers/*.h*") #Добавляем в переменную hdrs все файлы с расширением .h
file(GLOB srcs "sources/*.cpp") #Добавляем в переменную srcs все файлы с расширением .cpp

add_library(${alg} STATIC ${srcs} ${hdrs})