#ifndef __ALGORITHMS_HPP__
#define __ALGORITHMS_HPP__

#include <cmath>

#include "Field.hpp"

constexpr double PI = 3.14159265358979323846;
constexpr double C = 29979245800.0;  // speed of light (CGS)


namespace FDTD {

  class FDTD {
  public:
    FDTD() = delete;
    FDTD(const std::pair<uint64_t, uint64_t>& Nx_Ny,
      const std::pair<double, double>& ax_ay,
      const std::pair<double, double>& bx_by, double _dt);
    FDTD(const FDTD& _fields);
    ~FDTD() = default;


    uint64_t get_Nx(void) const;
    uint64_t get_Ny(void) const;

    Field::ComputingField& get_Ex(void) { return Ex; }
    Field::ComputingField& get_Ey(void) { return Ey; }
    Field::ComputingField& get_Ez(void) { return Ez; }
    Field::ComputingField& get_Bx(void) { return Bx; }
    Field::ComputingField& get_By(void) { return By; }
    Field::ComputingField& get_Bz(void) { return Bz; }

    void field_update(const double t);

  private:
    uint64_t Nx, Ny;
    Field::ComputingField Ex, Ey, Ez, Bx, By, Bz;
    double ax, bx, ay, by, dx, dy, dt;
  };

}  // namespace FDTD

#endif  // !__ALGORITHMS_HPP__
