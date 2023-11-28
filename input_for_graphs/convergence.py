import csv
import matplotlib.pyplot as plt
import numpy as np

with open("convergence.csv", "r") as file:
    reader = csv.reader(file, delimiter=" ")
    convergence = [np.float64(x) for x in next(reader)]

x_arr = []
for i in range(0, np.size(convergence)):
    x_arr.append(pow(2, i))

plt.title('График сходимости')
# plt.plot(x, my_data_Bx, label='Численное решение')
plt.plot(x, x_arr, label='Сходимость')
plt.xlabel('Множитель')
# plt.ylabel('Поле Bx')
plt.legend()
plt.grid()

plt.show()