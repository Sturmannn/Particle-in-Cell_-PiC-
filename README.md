# Particle-in-Cell_-PiC-

### Запуск модуля ***Tests*** (руководство пользователя) (УСТАРЕЛО):
* В первых тестах файла *tests.hpp* выполняются вычисления полей аналитическим и численным способами. После выполнения данных тестов происходит запись результатов в файлы *my_data.csv* и *analytical_data.csv*, находящихся в директории **input_for_graphs**.
* При запуске исполняемого файла *Test.exe* в консоли отобразятся порядки сходимости центральных схем со сдвигами и без сдвигов, с заданными начальными условиями, определёнными в исходном файле *tests.hpp*. Соответственно, меняя начальные данные, можно проверить результат сходимости схем.
* Для построения графиков необходимо запустить *main.py* после выполнения *Test.exe*.
