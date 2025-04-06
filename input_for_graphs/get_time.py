import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

path_to_measurements_file = Path(__file__).resolve().parents[1].joinpath("measurements.csv")
python_measurements_path = path_to_measurements_file.parent.joinpath("py_measurements.csv")

# Создание файла, если он не существует
path_to_measurements_file.touch(exist_ok=True)

print(path_to_measurements_file)

df = pd.read_csv(path_to_measurements_file, header=None, sep='\t')
df.columns = [
    'time', 
    'threads', 
    'Procs', 
    'Proc OХ', 
    'Proc OУ', 
    'Proc OZ'
]
num_runs = 5

print(df)

# Сортировка по количеству процессов по осям X, Y и Z
sorted_df = df.sort_values(by=['Proc OХ', 'Proc OУ', 'Proc OZ'])
print(sorted_df)

time_column = sorted_df['time']

sublist = [] # Для запуска 'num_runs' процессов
for i in range(0, len(time_column), num_runs):
    sublist.append(time_column.iloc[i:i+num_runs].to_list())


averages = []
# Удаление максимального и минимального значения из каждого подмассива
for row in sublist:
    if len(row) > 2:
        print(np.argmax(row))
        print(np.argmin(row))
        row = np.delete(row, np.argmax(row))
        row = np.delete(row, np.argmin(row))
        averages.append(np.mean(row))
    else:
        print("Error: Need more data!")

print(sublist, "averages")
print(averages)

# proc_list = [1, 2, 4, 6, 8, 16, 32]
# thread_list = [32, 16, 8, 6, 4, 2, 1]

proc_list = [1, 2, 4, 8, 16]
# thread_list = [16, 8, 4, 2, 1]

# Здесь происходит автоматическая перезапись файла
df_write = pd.DataFrame(averages)
df_write.to_csv(python_measurements_path, index=False, header=False)
print(df_write)

# Построение графиков
plt.figure(figsize=(10, 6))

time_16 = [0.81, 0.43, 0.11, 0.09, 0.1]
time_32 = [3.68, 1.91, 0.49, 0.46, 0.52]
time_64 = [25.26, 13.13, 3.14, 2.26, 2.50]
time_128 = [307.18, 158.82, 34.38, 21.63, 23.35]

new_time_16 = [0.81, 0.43, 0.11, 0.09, 0.1]
new_time_32 = [3.68, 1.91, 0.49, 0.46, 0.52]
new_time_64 = [25.26, 13.13, 3.14, 2.26, 2.50]
new_time_128 = [307.18, 158.82, 34.38, 21.63, 23.35]

log_time_16 = np.log(time_16)
log_time_32 = np.log(time_32)
log_time_64 = np.log(time_64)
log_time_128 = np.log(time_128)
# plt.plot(proc_list, averages, label="16x16x16")
# plt.plot(proc_list, log_time_16, label="16x16x16")
# plt.plot(proc_list, log_time_32, label="32x32x32")
# plt.plot(proc_list, log_time_64, label="64x64x64")
# plt.plot(proc_list, log_time_128, label="128x128x128")
# plt.ylabel('Логарифмическое время (сек)')

# Запуск 128х128х128 с одним потоком и 1,2,4,8,16 процессами на одном узле (Сильная масштабируемость)
time_128 = [209.23, 111.28, 20.97, 11.41, 11.91]
log_time_128 = np.log(time_128)
plt.plot(proc_list, log_time_128, marker='o', label='Сильная масштабируемость')
# plt.xlabel('Число процессов')
# plt.ylabel('Логарифмическое время (сек)')
# plt.title('128х128х128 без openMP (СИЛЬНАЯ МАСШТАБИРУЕМОСТЬ)')
# Слабая масштабируемость
time_weak_scale = [0.21, 0.21, 0.24, 0.22, 0.36]
log_time_weak_scale = np.log(time_weak_scale)
plt.plot(proc_list, time_weak_scale, marker='o', label='Слабая масштабируемость')
plt.title('Сильная и слабая масштабируемости (Логарифм)')
plt.xlabel('Число процессов')
plt.ylabel('Время log (сек)')

# Настройка осей и легенды
# plt.xlabel('Число процессов')
# # plt.xticks(proc_list)
# plt.title('MPI + OpenMP (16)')
plt.legend()

# Отображение графика
plt.grid(True)
plt.show()