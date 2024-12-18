#ifndef __FIELD_HPP__
#define __FIELD_HPP__

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <mpi.h>
#include <filesystem>

// constexpr const char *path_to_analytic_data = PATH_TO_ANALYTICAL_DATA;
// constexpr const char *path_to_calculated_data = PATH_TO_CALCULATED_DATA;
// constexpr const char *path_to_convergence_data = PATH_TO_CONVERGENCE_DATA;

constexpr const char *path_to_analytical_data_directory = PATH_TO_ANALYTICAL_DATA_DIRECTORY; 
constexpr const char *path_to_calculated_data_directory = PATH_TO_CALCULATED_DATA_DIRECTORY;
constexpr const char *path_to_convergence_data = PATH_TO_CONVERGENCE_DATA;

enum class Axis : int { Ox, Oy, Oz };

namespace fs = std::filesystem;

namespace Field {

class ComputingField {
public:
  // Конструкторы и деструктор
  ComputingField() = delete;
  ComputingField(const int64_t _Nx, const int64_t _Ny, const int64_t _Nz = 1);
  ComputingField(const ComputingField &_field);
  ComputingField(ComputingField &&_field) noexcept;
  ~ComputingField() = default;

  // Методы для получения размеров сетки
  int64_t get_Nx() const noexcept { return Nx; }
  int64_t get_Ny() const noexcept { return Ny; }
  int64_t get_Nz() const noexcept { return Nz; }

  std::vector<double> &get_field() { return field; }
  const std::vector<double> &get_field() const { return field; }

  double *data() { return field.data(); }
  int64_t fullsize() const { return field.size(); }

  void resize_field(const int64_t _Nx, const int64_t _Ny,
                    const int64_t _Nz = 1);

  // Операторы доступа к элементам поля
  double &operator()(const int64_t i, const int64_t j, const int64_t k = -1);
  const double &operator()(const int64_t i, const int64_t j,
                           const int64_t k = -1) const;

  // Операторы присваивания
  ComputingField &operator=(const ComputingField &_field);
  ComputingField &operator=(ComputingField &&_field) noexcept;

  // Запись поля в файл по фиксированной координате
  void write_field_to_file_OX(const char *path, const int64_t j = 0);
  void write_field_to_file_OY(const char *path, const int64_t i = 0);
  void write_field_to_file_OZ(const char *path, const int64_t j = 0);

  void clear_field() noexcept;

  // Очистка файла
  // static void clear_file(const char *path);

  // Очистка файлов
  static void clear_files(const char *directory_path);

private:
  void write_field_to_file(const char *path, const int64_t index, Axis axis);
  // Приватные члены данных
  std::vector<double> field; // Поле в виде одномерного вектора
  int64_t Nx, Ny, Nz; // Количество ячеек в сетке
};

} // namespace Field
#endif // !__FIELD_HPP__