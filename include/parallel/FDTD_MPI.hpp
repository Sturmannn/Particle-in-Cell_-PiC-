#pragma once

#include "../сore/types.hpp"
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <omp.h>
#include <string>

#include "../сore/Grid.hpp"
#include <memory>

namespace FDTD {

class FDTD_MPI {
public:
  FDTD_MPI(std::shared_ptr<Grid> _grid);
  FDTD_MPI(const FDTD_MPI &_fields) = default; // Copy constructor!
  ~FDTD_MPI() = default;

  const Field &EX(void) const { return Ex; }
  const Field &EY(void) const { return Ey; }
  const Field &EZ(void) const { return Ez; }
  const Field &BX(void) const { return Bx; }
  const Field &BY(void) const { return By; }
  const Field &BZ(void) const { return Bz; }

  Field &EX(void) { return Ex; }
  Field &EY(void) { return Ey; }
  Field &EZ(void) { return Ez; }
  Field &BX(void) { return Bx; }
  Field &BY(void) { return By; }
  Field &BZ(void) { return Bz; }

  FDTD_MPI &operator=(const FDTD_MPI &_fields) = default;

  std::shared_ptr<Grid> get_grid() const { return grid; }
  void Courant_condition_check(const Shift _shift) const;

  void shifted_field_update(const int64_t t, MPI_Comm cart_comm);

  static Axis get_axis(const Component E, const Component B);
  static std::string axisToString(const Component E, const Component B);

  void boundary_synchronization_3D();
  virtual void solve(const int t) = 0; // Pure virtual function

protected:
  std::shared_ptr<Grid> grid;

  MPI_Comm cart_comm;
  Field Ex, Ey, Ez, Bx, By, Bz;
  void update_E_field(const int iterations);
  void update_B_field(const int iterations);
};

class AnalyticalSolverFDTD : public FDTD_MPI {
public:
  using FDTD_MPI::FDTD_MPI; // Inherit constructor

  void solve(const int t) override {}

protected:
  void analytical_soulution(const Component E, const Component B,
                            const double t, const Shift _shift);
};

class NumericalSolverFDTD : public FDTD_MPI {
public:
  using FDTD_MPI::FDTD_MPI; // Inherit constructor

  void solve(const int t) override {}

protected:
};

} // namespace FDTD