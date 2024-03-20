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
  ComputingField() = delete;
  ComputingField(const uint64_t _Nx, const uint64_t _Ny, const uint64_t _Nz = 1ull);
  ComputingField(const ComputingField &_field);
  ComputingField(ComputingField &&_field) noexcept;
  ~ComputingField() = default;

  // Методы для получения размеров сетки
  uint64_t get_Nx() const noexcept { return Nx; }
  uint64_t get_Ny() const noexcept { return Ny; }
  uint64_t get_Nz() const noexcept { return Nz; }

  double* data() { return field.data(); }
  
  void resize_field(const uint64_t _Nx, const uint64_t _Ny, const uint64_t _Nz = 1ull);

  // Операторы доступа к элементам поля
  double &operator()(int64_t i, int64_t j, int64_t k = 0ull);
  const double &operator()(int64_t i, int64_t j, int64_t k = 0ull) const;

  // Операторы присваивания
  ComputingField &operator=(const ComputingField &_field);
  ComputingField &operator=(ComputingField &&_field) noexcept;

  // Запись поля в файл по фиксированной координате OX
  void write_field_to_file_OX(const char *path, const uint64_t j = 0ull);

  // Запись поля в файл по фиксированной координате OY
  void write_field_to_file_OY(const char* path, const uint64_t i = 0ull);

  void write_field_to_file_OZ(const char* path, uint64_t k = 0ull);

  // Очистка файла
  static void clear_file(const char *path);

  std::vector<double> field; // Поле в виде одномерного вектора
private:
  // Приватные члены данных
  uint64_t Nx, Ny, Nz;           // Количество ячеек в сетке
};

} // namespace Field
#endif // !__FIELD_HPP__