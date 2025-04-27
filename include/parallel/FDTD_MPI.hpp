#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <mpi.h>
#include <omp.h>
#include <string>

#include "Grid.hpp"
#include "MPI_Wrapper.hpp"
#include "types.hpp"

namespace FDTD {

class FDTD_MPI {
public:
  FDTD_MPI(std::shared_ptr<Grid> _grid,
           std::shared_ptr<MPI_Wrapper> _mpi_wrapper);
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

  double get_sign(const Component E, const Component B) const;
  Field &get_E_field(const Component E);
  Field &get_B_field(const Component B);

  void boundary_synchronization_3D();

  void shifted_field_update(const int64_t t, MPI_Comm cart_comm);

  static Axis get_axis(const Component E, const Component B);
  static std::string axisToString(const Component E, const Component B);

  double get_space_delta(const Component E, const Component B) const;

  virtual void solve(const Component E, const Component B, const double t,
                     const Shift _shift) = 0; // Pure virtual function

protected:
  std::shared_ptr<Grid> grid;
  std::shared_ptr<MPI_Wrapper> mpi_wrapper;

  MPI_Comm cart_comm;
  Field Ex, Ey, Ez, Bx, By, Bz;
  void update_E_field(const int iterations);
  void update_B_field(const int iterations);
};

class AnalyticalSolverFDTD : public FDTD_MPI {
public:
  using FDTD_MPI::FDTD_MPI; // Inherit constructor

  void solve(const Component E, const Component B, const double t,
             const Shift _shift) override;

protected:
  void analytical_soulution(const Component E, const Component B,
                            const double t, const Shift _shift);
};

class NumericalSolverFDTD : public FDTD_MPI {
public:
  using FDTD_MPI::FDTD_MPI; // Inherit constructor
                            // double t
  void solve(const Component E, const Component B, const double t,
             const Shift _shift) override;

protected:
  void update_E_field();
  void update_B_field();
  void set_default_values(const Component E, const Component B,
                          const Shift _shift);
  void numerical_solution(const double t);
};

} // namespace FDTD