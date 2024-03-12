#ifndef __FDTD_HPP__
#define __FDTD_HPP__

#include <cmath>
#include <tuple>

#include "Field.hpp"

constexpr double PI = 3.14159265358979323846;
constexpr double C = 29979245800.0; // speed of light (CGS)

namespace FDTD {

class FDTD {
public:
  // Конструкторы и деструктор
  FDTD() = delete;
  FDTD(const std::tuple<uint64_t, uint64_t, uint64_t> &Nx_Ny_Nz,
       const std::tuple<double, double, double> &ax_ay_az,
       const std::tuple<double, double, double> &bx_by_bz, double _dt);
  FDTD(const FDTD &_fields);
  FDTD(FDTD &&_fields) noexcept;
  ~FDTD() = default;

  // Операторы присваивания
  FDTD &operator=(const FDTD &_fields);
  FDTD &operator=(FDTD &&_fields) noexcept;

  // Методы для получения размеров сетки и коэффициентов
  uint64_t get_Nx(void) const noexcept { return Nx; }
  uint64_t get_Ny(void) const noexcept { return Ny; }
  uint64_t get_Nz(void) const noexcept { return Nz; }

  std::pair<double, double> get_ax_bx(void) const noexcept {
    return std::make_pair(ax, bx);
  }
  std::pair<double, double> get_ay_by(void) const noexcept {
    return std::make_pair(ay, by);
  }
  std::pair<double, double> get_az_bz(void) const noexcept {
    return std::make_pair(az, bz);
  }

  // Методы для получения компонент полей
  Field::ComputingField &get_Ex(void) { return Ex; }
  Field::ComputingField &get_Ey(void) { return Ey; }
  Field::ComputingField &get_Ez(void) { return Ez; }
  Field::ComputingField &get_Bx(void) { return Bx; }
  Field::ComputingField &get_By(void) { return By; }
  Field::ComputingField &get_Bz(void) { return Bz; }

  // Методы для получения параметров сетки
  double get_dx(void) const noexcept { return dx; }
  double get_dy(void) const noexcept { return dy; }
  double get_dz(void) const noexcept { return dz; }
  double get_dt(void) const noexcept { return dt; }

  void set_dt(double _dt) { dt = _dt; }
  void set_Nx_Ny_Nz(uint64_t _Nx, uint64_t _Ny, uint64_t _Nz);

  // Методы для обновления полей
  void field_update(const double t);
  void field_update(const uint64_t t);

  void shifted_field_update(const double t);
  void shifted_field_update(const uint64_t t);

  // Запись полей в файлы по фиксированным координатам OX и OY
  void write_fields_to_file_OX(const char *path, const double dx,
                               uint64_t j = 0ull); // The row is fixed
  void write_fields_to_file_OY(const char *path, const double dy,
                               uint64_t i = 0ull); // The col is fixed
  void write_fields_to_file_OZ(const char* path, const double dz,
    uint64_t k = 0ull); // The col is fixed
private:
  // Приватные члены данных
  uint64_t Nx, Ny, Nz;                          // Размеры сетки
  Field::ComputingField Ex, Ey, Ez, Bx, By, Bz; // Компоненты полей
  double ax, bx, ay, by, az, bz, dx, dy, dz, dt; // Коэффициенты и шаги
  uint64_t check_Nz_dimension(uint64_t _Nz) { return (_Nz > 1 ?  _Nz :  0); }
};

} // namespace FDTD

#endif // !__FDTD_HPP__
