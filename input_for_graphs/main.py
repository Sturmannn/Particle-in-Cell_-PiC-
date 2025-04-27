import os
import numpy as np
import matplotlib.pyplot as plt
import glob

def check_convergence():
    print("Check convergence")

plt.style.use('ggplot')

# Получаем путь к директории, где находится текущий файл
file_directory = os.path.dirname(os.path.abspath(__file__))

# Объединение текущей директории с указанными путями
numerical_path  = os.path.join(file_directory, 'my_data')
analytical_path  = os.path.join(file_directory, 'analytical_data')

fields = ['Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz']

numerical_data = {field: [] for field in fields}
analytical_data = {field: [] for field in fields}
deltas = []

# Считывание всех файлов с данными
numerical_files = sorted(glob.glob(os.path.join(numerical_path, 'my_data_*.csv')))
analytical_files = sorted(glob.glob(os.path.join(analytical_path, 'analytical_data_*.csv')))

num_combined_data = []
anl_combined_data = []
delta = None

# Чтение и объединение данных численного и аналитического решений
for num_file, anl_file in zip(numerical_files, analytical_files):

    num_data = np.genfromtxt(num_file, delimiter=';', skip_footer=1)
    anl_data = np.genfromtxt(anl_file, delimiter=';', skip_footer=1)

    # Добавление данных к существующим строкам
    if not num_combined_data:
        num_combined_data = num_data.tolist()
        anl_combined_data = anl_data.tolist()
    else:
        for i in range(len(num_combined_data)):
            num_combined_data[i].extend(num_data[i])
            anl_combined_data[i].extend(anl_data[i])
    
    # Чтение дельты, если она еще не была считана
    if delta is None:
        delta = np.genfromtxt(num_file, delimiter=';', skip_header=len(num_data))
        if isinstance(delta, np.ndarray):
            delta = delta.item()  # Преобразование массива в число

num_combined_data = np.array(num_combined_data, dtype=object)
anl_combined_data = np.array(anl_combined_data, dtype=object)

# Создание графиков для каждого поля
fig, axs = plt.subplots(3, 2, figsize=(14, 10))
axs = axs.flatten()

# Создание оси Ox
x_axis = np.arange(len(anl_combined_data[0]))
# x_axis = np.arange(len(num_combined_data[0])) * delta

# Цвета для различных графиков
colors = ['blue', 'red']

# Построение графиков
for idx, field in enumerate(fields):
    # Получение данных для численного и аналитического решения
    num_values = num_combined_data[idx]
    anl_values = anl_combined_data[idx]
    
    # Построение численного и аналитического графиков
    axs[idx].plot(x_axis, anl_values, label='Аналитическое решение', color=colors[1], linewidth=3, alpha=0.7)
    # axs[idx].plot(x_axis, num_values, label='Численное решение', color=colors[0], linestyle='--', linewidth=4, alpha=0.7)
    
    # Настройка подписей и заголовков
    axs[idx].set_xlabel(r'N$\Delta$x', fontsize=12)
    axs[idx].set_ylabel(f"Поле {field}", fontsize=12)
    axs[idx].set_title(f"График значения {field}", fontsize=14)
    
    # Добавление легенды и сетки
    axs[idx].legend(fontsize=10)
    axs[idx].grid(True, color='black', linestyle=':')
    
    # # Увеличиваем шрифт
    # for label in (axs[idx].get_xticklabels() + axs[idx].get_yticklabels()):
    #     label.set_fontsize(10)

plt.tight_layout()
plt.show()
