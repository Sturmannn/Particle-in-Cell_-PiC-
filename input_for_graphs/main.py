import csv
import matplotlib.pyplot as plt
import numpy as np

with open("analytical_data.csv", "r") as file:
    reader = csv.reader(file, delimiter=";")
    an_data_first_row = next(reader)
    an_data = [np.float64(x) for x in an_data_first_row]
    an_data_second_row = next(reader)
    an_data_dx = np.float64(an_data_second_row[0])

with open("my_data.csv", "r") as file:
    reader = csv.reader(file, delimiter=";")
    my_data_first_row = next(reader)
    my_data = [np.float64(x) for x in my_data_first_row]
    my_data_second_row = next(reader)
    my_data_dx = np.float64(my_data_second_row[0])

dx = 0
x = []
for i in range(0, np.size(my_data)):
    x.append(i + dx)
    dx += my_data_dx

plt.xlabel("dx")
plt.ylabel("Значение Еу")
plt.title('Графики значения Ey')

plt.plot(x, an_data, label="Аналитический график")
plt.plot(x, my_data, label="Численный график")


plt.legend()
plt.grid()
plt.show()
