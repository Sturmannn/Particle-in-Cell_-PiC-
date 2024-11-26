file(GLOB hdrs "headers/*.h*" "*.hpp") #Добавляем в переменную hdrs все файлы с расширением .h .hpp
file(GLOB srcs "sources/*.cpp" "sources/*.cc") #Добавляем в переменную srcs все файлы с расширением .cpp .cc

add_library(${TEST_LIB} STATIC ${srcs} ${hdrs}) #Создаем статическую библиотеку с именем из переменной 