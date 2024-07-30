#ifndef __FDTD_HPP__
#define __FDTD_HPP__

#include <cmath>
#include <mpi.h>
#include <tuple>

#include "Field.hpp"

//#define MPI

constexpr double PI = 3.14159265358979323846;
constexpr double C = 29979245800.0; // speed of light (CGS)

namespace FDTD {

enum class Component { Ex, Ey, Ez, Bx, By, Bz };
enum class Axis { Ox, Oy, Oz };

class FDTD {
public:
  // Конструкторы и деструктор
  FDTD() = delete;
  FDTD(const std::tuple<int64_t, int64_t, int64_t> &Nx_Ny_Nz,
       const std::tuple<double, double, double> &ax_ay_az,
       const std::tuple<double, double, double> &bx_by_bz, double _dt);
  FDTD(const FDTD &_fields);
  FDTD(FDTD &&_fields) noexcept;
  ~FDTD() = default;

  // Операторы присваивания
  FDTD &operator=(const FDTD &_fields);
  FDTD &operator=(FDTD &&_fields) noexcept;

  // Методы для получения размеров сетки и коэффициентов
  int64_t get_Nx(void) const noexcept { return Nx; }
  int64_t get_Ny(void) const noexcept { return Ny; }
  int64_t get_Nz(void) const noexcept { return Nz; }

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
  void set_Nx_Ny_Nz(int64_t _Nx, int64_t _Ny, int64_t _Nz);

  // Методы для обновления полей
  void field_update(const double t);
  void field_update(const int64_t t);

  void shifted_field_update(const double t);
  void shifted_field_update(const int64_t t);

  void
  write_fields_to_file(const char *path, Component E, Component B,
                       const double delta,
                       const int64_t row_number = 0ull); // The col is fixed

  static Axis get_axis(const Component E, const Component B);

private:
  // Приватные члены данных
  int64_t Nx, Ny, Nz;                           // Размеры сетки
  Field::ComputingField Ex, Ey, Ez, Bx, By, Bz; // Компоненты полей
  double ax, bx, ay, by, az, bz, dx, dy, dz, dt; // Коэффициенты и шаги
};

} // namespace FDTD

#endif // !__FDTD_HPP__
