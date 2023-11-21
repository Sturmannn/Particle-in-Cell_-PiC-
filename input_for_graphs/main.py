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

with open("analytical_data.csv", "r") as file:
    reader = csv.reader(file, delimiter=";")
    an_data_Ex = [np.float64(x) for x in next(reader)]
    an_data_Ey = [np.float64(x) for x in next(reader)]
    an_data_Ez = [np.float64(x) for x in next(reader)]
    an_data_Bx = [np.float64(x) for x in next(reader)]
    an_data_By = [np.float64(x) for x in next(reader)]
    an_data_Bz = [np.float64(x) for x in next(reader)]
    an_data_dx = np.float64(next(reader)[0])

with open("my_data.csv", "r") as file:
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

# plt.xlabel('$\Delta x$')
# plt.ylabel("Поле Еу")
# plt.title('Графики значения Ey')

# print(len(my_data_Ey))
# print(np.size(my_data_dx))

# plt.plot(x, an_data_Ey, label="Аналитический график")
# plt.plot(x, my_data_Ey, label="Численный график")
# # plt.savefig(".\Plots\myPlot.png")

# # my_div_an = np.asarray(my_data) / np.asarray(an_data);
# # plt.plot(x, my_div_an, label="График деления АН на ЧИ");

# plt.legend()
# plt.grid()
# plt.show()

# Создание шести осей в трех рядах и два столбца
fig, axs = plt.subplots(3, 2)

# Построение графиков на каждой из созданных осей
# Поля E
axs[0, 0].plot(x, my_data_Ex, label='Численное решение')
axs[0, 0].plot(x, an_data_Ex, label='Аналитическое решение')
axs[0, 0].set_xlabel(r'N$\Delta$x')
axs[0, 0].set_ylabel('Поле Ex')
axs[0, 0].legend()
axs[0, 0].grid()

axs[1, 0].plot(x, my_data_Ey, label='Численное решение')
axs[1, 0].plot(x, an_data_Ey, label='Аналитическое решение')
axs[1, 0].set_xlabel(r'N$\Delta$x')
axs[1, 0].set_ylabel('Поле Ey')
axs[1, 0].legend()
axs[1, 0].grid()

axs[2, 0].plot(x, my_data_Ez, label='Численное решение')
axs[2, 0].plot(x, an_data_Ez, label='Аналитическое решение')
axs[2, 0].set_xlabel(r'N$\Delta$x')
axs[2, 0].set_ylabel('Поле Ez')
axs[2, 0].legend()
axs[2, 0].grid()

# Поля B
axs[0, 1].plot(x, my_data_Bx, label='Численное решение')
axs[0, 1].plot(x, an_data_Bx, label='Аналитическое решение')
axs[0, 1].set_xlabel(r'N$\Delta$x')
axs[0, 1].set_ylabel('Поле Bx')
axs[0, 1].legend()
axs[0, 1].grid()

axs[1, 1].plot(x, my_data_By, label='Численное решение')
axs[1, 1].plot(x, an_data_By, label='Аналитическое решение')
axs[1, 1].set_xlabel(r'N$\Delta$x')
axs[1, 1].set_ylabel('Поле By')
axs[1, 1].legend()
axs[1, 1].grid()

axs[2, 1].plot(x, my_data_Bz, label='Численное решение')
axs[2, 1].plot(x, an_data_Bz, label='Аналитическое решение')
axs[2, 1].set_xlabel(r'N$\Delta$x')
axs[2, 1].set_ylabel('Поле Bz')
axs[2, 1].legend()
axs[2, 1].grid()

# plt.plot(x, my_data_By, label='Численное решение')
# plt.plot(x, an_data_By, label='Аналитическое решение')
# plt.xlabel(r'N$\Delta$x')
# plt.ylabel('Поле By')
# plt.legend()
# plt.grid()


plt.show()
