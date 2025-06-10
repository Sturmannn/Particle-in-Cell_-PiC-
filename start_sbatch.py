# -*- coding: utf-8 -*-

import subprocess
import time
import os

from pathlib import Path

# path_to_measurements_file = Path(__file__).resolve().parent.joinpath("measurements.csv")

# python_measurements_path = path_to_measurements_file.parent.joinpath("py_measurements.csv")

current_file_path = os.path.abspath(__file__)
current_dir = os.path.dirname(current_file_path)

path_to_measurements_file = os.path.join(current_dir, "measurements.csv")
python_measurements_path = os.path.join(os.path.dirname(path_to_measurements_file), "py_measurements.csv")


# Очистка файла с измерениями
with open(path_to_measurements_file, "w"):
    pass

# path_to_bin_file = Path(__file__).resolve().parent.joinpath("bin", "Tests")
# path_to_bin_file = current_dir.join("bin", "Tests")
path_to_bin_file = os.path.join(current_dir, "bin", "SampleFDTD")

# path_to_bin_file.joinpath("bin", "Tests")
print('Путь до исполняемого файла -', path_to_bin_file)

num_threads = 1
num_runs = 5 # Количество итераций запуска программы для подсчёта времени работы
ntasks_per_node = 1 # Число процессов на узел
num_clusters = 1 # Число узлов
np = num_clusters * num_clusters # Общее число процессов
t = 200 # Лимит на выполнение задачи

# while True:
#     try:
#         user_input = input('Введите число запусков приложения: ')
#         num_runs = int(user_input)

#         user_input = input('Введите число процессов на узел: ')
#         ntasks_per_node = int(user_input)

#         user_input = input('Введите общее число узлов: ')
#         num_clusters = int(user_input)

#         np = ntasks_per_node * num_clusters
#         break
#     except ValueError:
#         print('Ошибка! Вводите корректные числа')

# def call_sbatch(_time, _ntasks_per_node, _num_clusters, _path_to_bin_file, _OMP_NUM_THREADS=1):
#     sbatch_command = 'sbatch -v -t {_time} -p gpu --ntasks-per-node {_ntasks_per_node} -N {_num_clusters} \
#           --wrap="export OMP_NUM_THREADS={_OMP_NUM_THREADS}; mpiexec -np {_ntasks_per_node * _num_clusters} {_path_to_bin_file}"' \
#           .format(_time, _ntasks_per_node, _num_clusters, _OMP_NUM_THREADS, _ntasks_per_node * _num_clusters, _path_to_bin_file)
#     return sbatch_command

# def call_sbatch(_time, _ntasks_per_node, _num_clusters, _path_to_bin_file, _domain_size,  _OMP_NUM_THREADS=1):
#     sbatch_command = 'sbatch -v -t {} -p gpu --ntasks-per-node {} -N {} \
#           --wrap="export OMP_NUM_THREADS={}; mpiexec -n {} {} --Nx {} --Ny {} --Nz {}"' \
#           .format(_time, _ntasks_per_node, _num_clusters, \
#                     _OMP_NUM_THREADS, _ntasks_per_node * _num_clusters, \
#                     _path_to_bin_file, _domain_size, _domain_size, _domain_size)
#     return sbatch_command

def call_sbatch(_time, _ntasks_per_node, _num_clusters, _path_to_bin_file, _domain_size, _OMP_NUM_THREADS=1):
    total_tasks = _ntasks_per_node * _num_clusters
    sbatch_command = (
        f'sbatch -v -t {_time} -p gpu --ntasks-per-node {_ntasks_per_node} -N {_num_clusters} '
        f'--wrap="export OMP_NUM_THREADS={_OMP_NUM_THREADS}; '
        f'mpiexec -n {total_tasks} {_path_to_bin_file}"'
    )
    return sbatch_command



# Устанавливаем количество OpenMP потоков
# os.environ["OMP_NUM_THREADS"] = num_threads


# print(f'sbatch -v -t {t} -p gpu --ntasks-per-node {ntasks_per_node} -N {num_clusters} --wrap="mpiexec -np {np} {path_to_bin_file}"')


domain_sizes = [512]
# thread_list = [1, 2, 4, 8, 16, 32]
# thread_list = [1, 1, 1, 1, 1, 1, 1, 1]
# thread_list = [1, 1, 1, 1, 1, 1]
thread_list = [1, 2, 4, 8, 16, 32]
proc_list = [1, 1, 1, 1, 1, 1]
# proc_list = [1, 2, 4, 8, 16, 32]
# proc_list = [32, 16, 8, 4, 2, 1]
# proc_list = [1, 2, 3, 4, 5, 6, 7, 8]


# for index in range(len(thread_list)):
for sizes in range(len(domain_sizes)):
    for count in range(num_runs):
        for i in range(len(proc_list)):
            print('Запуск #{}'.format(i))
            try:
                # subprocess.run(call_sbatch(t, proc_list[index], num_clusters, path_to_bin_file, thread_list[index]), \
                #                 shell=True, check=True)
                subprocess.run(call_sbatch(t, proc_list[i], num_clusters, path_to_bin_file, domain_sizes[sizes], thread_list[i]), \
                        shell=True, check=True)
                time.sleep(2)
            except subprocess.SubprocessError as er:
                print('Запуск скрипта выполнен с ошибкой: {}'.format(er))
        


# start_time = time.time()
# end_time = time.time()


# elapsed_time = end_time - start_time
# print(f'Время работы скрипта - {elapsed_time} сек.')