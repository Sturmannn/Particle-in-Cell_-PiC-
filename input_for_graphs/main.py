import csv
import matplotlib.pyplot as plt
import numpy as np

# # with open("analytical_data.csv", "r") as file:
# #     reader = csv.reader(file, delimiter=";")
# #     an_data_first_row = next(reader)
# #     an_data = [np.float64(x) for x in an_data_first_row]
# #     an_data_second_row = next(reader)
# #     an_data_dx = np.float64(an_data_second_row[0])

# # with open("my_data.csv", "r") as file:
# #     reader = csv.reader(file, delimiter=";")
# #     my_data_first_row = next(reader)
# #     my_data = [np.float64(x) for x in my_data_first_row]
# #     my_data_second_row = next(reader)
# #     my_data_dx = np.float64(my_data_second_row[0])

with open("./input_for_graphs/analytical_data.csv", "r") as file:
    reader = csv.reader(file, delimiter=";")
    an_data_Ex = [np.float64(x) for x in next(reader)]
    an_data_Ey = [np.float64(x) for x in next(reader)]
    an_data_Ez = [np.float64(x) for x in next(reader)]
    an_data_Bx = [np.float64(x) for x in next(reader)]
    an_data_By = [np.float64(x) for x in next(reader)]
    an_data_Bz = [np.float64(x) for x in next(reader)]
    an_data_dx = np.float64(next(reader)[0])

with open("./input_for_graphs/my_data.csv", "r") as file:
    reader = csv.reader(file, delimiter=";")
    my_data_Ex = [np.float64(x) for x in next(reader)]
    my_data_Ey = [np.float64(x) for x in next(reader)]
    my_data_Ez = [np.float64(x) for x in next(reader)]
    my_data_Bx = [np.float64(x) for x in next(reader)]
    my_data_By = [np.float64(x) for x in next(reader)]
    my_data_Bz = [np.float64(x) for x in next(reader)]
    my_data_dx = np.float64(next(reader)[0])

dx = 0
x = []
for i in range(0, np.size(my_data_Ex)):
    x.append(i + dx)
    dx += my_data_dx

# plt.xlabel('N$\Delta x$', fontsize=20)
# plt.ylabel("Поле Еz", fontsize=20)
# plt.title('Графики значения Ez')
# plt.title('График значения Ez', fontsize=20)
# plt.title('Разность численного и аналитических решений (Ez)', fontsize=20)

print(len(my_data_Ez))
print(np.size(my_data_dx))

# sub_data = []
# for i in range(0, np.size(an_data_Ez)):
#     sub_data.append(my_data_Ez[i] - an_data_Ez[i])
# plt.plot(x, sub_data, linestyle='-', linewidth=3 , color='black', label="Разность решений", alpha=1)
# plt.plot(x, an_data_Ez, linestyle='-', linewidth=3 , color='r', label="Аналитическое решение", alpha=1)
# plt.plot(x, my_data_Ez, linestyle='--', linewidth=6, color='b', label="Численное решение", alpha=0.6)

# plt.plot(x, an_data_Ez, 'g', label="Аналитический график")
# plt.plot(x, my_data_Ez, label="Численный график")
# plt.savefig(".\Plots\myPlot.png")

# my_div_an = np.asarray(my_data) / np.asarray(an_data);
# plt.plot(x, my_div_an, label="График деления АН на ЧИ");

# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.legend(fontsize=20)
# plt.grid()
# plt.show()



# # Создание шести осей в трех рядах и два столбца
fig, axs = plt.subplots(3, 2)

# Построение графиков на каждой из созданных осей
# Поля E
axs[0, 0].plot(x, my_data_Ex, linestyle='--', linewidth=4, color="red", label='Численное решение')
axs[0, 0].plot(x, an_data_Ex, alpha=0.6, linewidth=3, label='Аналитическое решение')
axs[0, 0].set_xlabel(r'N$\Delta$x')
axs[0, 0].set_ylabel('Поле Ex')
axs[0, 0].legend()
axs[0, 0].grid()

axs[1, 0].plot(x, my_data_Ey, linestyle='--', linewidth=4, color="red", label='Численное решение')
axs[1, 0].plot(x, an_data_Ey, alpha=0.6, linewidth=3, label='Аналитическое решение')
axs[1, 0].set_xlabel(r'N$\Delta$x')
axs[1, 0].set_ylabel('Поле Ey')
axs[1, 0].legend()
axs[1, 0].grid()

axs[2, 0].plot(x, my_data_Ez, linestyle='--', linewidth=4, color="red", label='Численное решение')
axs[2, 0].plot(x, an_data_Ez, alpha=0.6, linewidth=3, label='Аналитическое решение')
axs[2, 0].set_xlabel(r'N$\Delta$x')
axs[2, 0].set_ylabel('Поле Ez')
axs[2, 0].legend()
axs[2, 0].grid()

# Поля B
axs[0, 1].plot(x, my_data_Bx, linestyle='--', linewidth=4, color="red", label='Численное решение')
axs[0, 1].plot(x, an_data_Bx, alpha=0.6, linewidth=3, label='Аналитическое решение')
axs[0, 1].set_xlabel(r'N$\Delta$x')
axs[0, 1].set_ylabel('Поле Bx')
axs[0, 1].legend()
axs[0, 1].grid()

axs[1, 1].plot(x, my_data_By, linestyle='--', linewidth=4, color="red", label='Численное решение')
axs[1, 1].plot(x, an_data_By, alpha=0.6, linewidth=3, label='Аналитическое решение')
axs[1, 1].set_xlabel(r'N$\Delta$x')
axs[1, 1].set_ylabel('Поле By')
axs[1, 1].legend()
axs[1, 1].grid()

axs[2, 1].plot(x, my_data_Bz, linestyle='--', linewidth=4, color="red", label='Численное решение')
axs[2, 1].plot(x, an_data_Bz, alpha=0.6, linewidth=3, label='Аналитическое решение')
axs[2, 1].set_xlabel(r'N$\Delta$x')
axs[2, 1].set_ylabel('Поле Bz')
axs[2, 1].legend()
axs[2, 1].grid()

# plt.plot(x, my_data_Bx, label='Численное решение')
# plt.plot(x, an_data_Bx, label='Аналитическое решение')
# plt.xlabel(r'N$\Delta$x')
# plt.ylabel('Поле Bx')
# plt.legend()
# plt.grid()


plt.show()
