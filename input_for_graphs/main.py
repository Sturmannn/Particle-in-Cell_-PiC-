import csv

name_1 = "Petya"
name_2 = "Gabov"

array_1 = [name_1, name_2]
array_2 = [name_1, name_1]

with open("data.csv", "w") as file:
    writer = csv.writer(file, delimiter=";")
    writer.writerows(
        [array_1, array_2, array_2]
    )