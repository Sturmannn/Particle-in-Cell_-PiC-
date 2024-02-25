#ifndef __FIELD_HPP__
#define __FIELD_HPP__

#include <fstream>
#include <iostream>
#include <omp.h>
#include <vector>


constexpr char *path_to_analytic_data =
    "..\\..\\input_for_graphs\\analytical_data.csv"; // Writing E, B to a file
constexpr char *path_to_calculated_data =
    "..\\..\\input_for_graphs\\my_data.csv";

constexpr char *path_to_convergence_data =
    "..\\..\\input_for_graphs\\convergence.csv"; // Writing convergence

namespace Field {

class ComputingField {
public:
  // Конструкторы и деструктор
  ComputingField(const uint64_t _Nx, const uint64_t _Ny);
  ComputingField(const ComputingField &_field);
  ComputingField(ComputingField &&_field) noexcept;
  ~ComputingField() = default;

  // Методы для получения размеров сетки
  uint64_t get_Nx() const noexcept { return Nx; }
  uint64_t get_Ny() const noexcept { return Ny; }

  // Операторы доступа к элементам поля
  double &operator()(uint64_t i, uint64_t j);
  const double &operator()(uint64_t i, uint64_t j) const;

  // Операторы присваивания
  ComputingField &operator=(const ComputingField &_field);
  ComputingField &operator=(ComputingField &&_field) noexcept;

  // Запись поля в файл по фиксированной координате OX
  void write_field_to_file_OX(const char *path, const uint64_t j = 0ull);

  // Запись поля в файл по фиксированной координате OY
  void write_field_to_file_OY(const char *path, const uint64_t i = 0ull);

  // Очистка файла
  static void clear_file(const char *path);

private:
  // Приватные члены данных
  uint64_t Nx, Ny;           // Количество ячеек в сетке
  std::vector<double> field; // Поле в виде одномерного вектора
};

} // namespace Field
#endif // !__FIELD_HPP__