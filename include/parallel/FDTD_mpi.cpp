#include "FDTD_MPI.hpp"

namespace FDTD {
inline FDTD_MPI::FDTD_MPI(std::shared_ptr<Grid> _grid,
                          std::shared_ptr<MPI_Wrapper> _mpi_wrapper)
    : grid(std::move(_grid)), mpi_wrapper(std::move(_mpi_wrapper)) {
  int grid_size = grid->mpi_get_total_bounded_subdomain_sizes();
  Ex = Ey = Ez = Bx = By = Bz = Field(grid_size);
}

inline void FDTD_MPI::Courant_condition_check(const Shift _shift) const {
  GridCoordinatesSteps steps = grid->get_coordinates_steps();
  TimeStep dt = grid->get_dt();
  double dx = steps.dx, dy = steps.dy, dz = steps.dz;
  if (_shift == Shift::shifted) {
    if (dt <= (dx / (C * sqrt(2))) && dt <= (dy / (C * sqrt(2))))
      return;
  } else {
    if (dt <= (dx * dx / 2 * (C * sqrt(2))) &&
        dt <= (dy * dy / 2 * (C * sqrt(2))))
      return;
  }
  std::cout << "Courant's condition is not satisfied\n";
  exit(-1);
}

inline int FDTD_MPI::get_sign(const Component E, const Component B) const {
  // Ox
  if (E == Component::Ey && B == Component::Bz)
    return 1;
  else if (E == Component::Ez && B == Component::By)
    return -1;
  // Oy
  else if (E == Component::Ez && B == Component::Bx)
    return 1;
  else if (E == Component::Ex && B == Component::Bz)
    return -1;
  // Oz
  else if (E == Component::Ex && B == Component::By)
    return 1;
  else if (E == Component::Ey && B == Component::Bx)
    return -1;
  else {
    std::cout << "\nGet_Sign: Error! Wrong Components! "
              << "E = " << E << " B = " << B << std::endl;
    exit(-1);
  }
}

inline Field &FDTD_MPI::get_E_field(const Component E) {
  if (E == Component::Ex)
    return Ex;
  else if (E == Component::Ey)
    return Ey;
  else if (E == Component::Ez)
    return Ez;
  else {
    std::cout << "\nGet_E_Field: Error! Wrong E - field!\n";
    exit(-1);
  }
}

inline Field &FDTD_MPI::get_B_field(const Component B) {
  if (B == Component::Bx)
    return Bx;
  else if (B == Component::By)
    return By;
  else if (B == Component::Bz)
    return Bz;
  else {
    std::cout << "\nGet_B_Field: Error! Wrong B - field!\n";
    exit(-1);
  }
}

inline Axis FDTD_MPI::get_axis(const Component E, const Component B) {
  if ((E == Component::Ey && B == Component::Bz) ||
      (E == Component::Ez && B == Component::By))
    return Axis::X;
  else if ((E == Component::Ez && B == Component::Bx) ||
           (E == Component::Ex && B == Component::Bz))
    return Axis::Y;
  else if ((E == Component::Ex && B == Component::By) ||
           (E == Component::Ey && B == Component::Bx))
    return Axis::Z;
  else {
    std::cout << "\nGet_Axis: Error! Wrong Components! "
              << "E = " << E << " B = " << B << std::endl;
    exit(-1);
  }
}

inline std::string FDTD_MPI::axisToString(const Component E,
                                          const Component B) {
  if (E == Component::Ey && B == Component::Bz)
    return "+ Ox";
  else if ((E == Component::Ez && B == Component::By))
    return "- Ox";
  else if ((E == Component::Ez && B == Component::Bx))
    return "+ Oy";
  else if ((E == Component::Ex && B == Component::Bz))
    return "- Oy";
  else if ((E == Component::Ex && B == Component::By))
    return "+ Oz";
  else if ((E == Component::Ey && B == Component::Bx))
    return "- Oz";
  else {
    std::cout << "\nGet_Axis: Error! Wrong Components! "
              << "E = " << E << " B = " << B << std::endl;
    exit(-1);
  }
}

void AnalyticalSolverFDTD::analytical_soulution(const Component E,
                                                const Component B,
                                                const double t,
                                                const Shift _shift) {
  Courant_condition_check(_shift);
  double coeff = (_shift == Shift::shifted) ? 0.5 : 0.0;

  Field &Ex = EX();
  Field &Ey = EY();
  Field &Ez = EZ();
  Field &Bx = BX();
  Field &By = BY();
  Field &Bz = BZ();

  Axis axis = get_axis(E, B);
  Field &E_field = get_E_field(E);
  Field &B_field = get_B_field(B);

  std::array<int, 3> neigbor_coords = {0, 0, 0}; // size = MPI_DIMENSION
  std::array<int, 3> coords = mpi_wrapper->get_coords();
  GridCoordinatesBounds bounds = grid->get_bounds();
  GridCoordinatesSteps steps = grid->get_coordinates_steps();
  int shift_relative_to_other_processes = 0;

  int Nx = grid->get_mpi_local_unbounded_subdomain_sizes().Nx;
  int Ny = grid->get_mpi_local_unbounded_subdomain_sizes().Ny;
  int Nz = grid->get_mpi_local_unbounded_subdomain_sizes().Nz;

  if (axis == Axis::X) {
    std::vector<int> all_local_Nx = mpi_wrapper->get_sizes_OX();
    for (int i = coords[0] - 1; coords[0] > 0; --i) {
      shift_relative_to_other_processes +=
          all_local_Nx[mpi_wrapper->get_cart_rank(i, coords[1], coords[2])];
    }
    int x = bounds.ax + steps.dx * shift_relative_to_other_processes;
    int sign = get_sign(E, B);
    int Nx = grid->get_mpi_local_unbounded_subdomain_sizes().Nx;

    for (int i = 0; i < Nx; ++i, x += steps.dx)
      for (int j = 0; j < Ny; ++j)
        for (int k = 0; k < Nz; ++k) {
          int index = Nx * Ny * (k + 1) + Nx * (j + 1) + (i + 1);
          Ex[index] = Ey[index] = Ez[index] = 0.0;
          Bx[index] = By[index] = Bz[index] = 0.0;

          E_field[index] = sin(2.0 * PI * (x - bounds.ax - sign * C * t) /
                               (bounds.bx - bounds.ax));
          B_field[index] =
              sin(2.0 * PI * (x + steps.dx * coeff - bounds.ax - sign * C * t) /
                  (bounds.bx - bounds.ax));
        }
  } else if (axis == Axis::Y) {
    std::vector<int> all_local_Ny = mpi_wrapper->get_sizes_OY();
    for (int i = coords[1] - 1; coords[1] > 0; --i) {
      shift_relative_to_other_processes +=
          all_local_Ny[mpi_wrapper->get_cart_rank(coords[0], i, coords[2])];
    }
    int y = bounds.ay + steps.dy * shift_relative_to_other_processes;
    int sign = get_sign(E, B);
    int Ny = grid->get_mpi_local_unbounded_subdomain_sizes().Ny;

    for (int i = 0; i < Nx; ++i, y += steps.dy)
      for (int j = 0; j < Ny; ++j)
        for (int k = 0; k < Nz; ++k) {
          int index = Nx * Ny * (k + 1) + Nx * (j + 1) + (i + 1);
          Ex[index] = Ey[index] = Ez[index] = 0.0;
          Bx[index] = By[index] = Bz[index] = 0.0;

          E_field[index] = sin(2.0 * PI * (y - bounds.ay - sign * C * t) /
                               (bounds.by - bounds.ay));
          B_field[index] =
              sin(2.0 * PI * (y + steps.dy * coeff - bounds.ay - sign * C * t) /
                  (bounds.by - bounds.ay));
        }
  } else if (axis == Axis::Z) {
    std::vector<int> all_local_Nz = mpi_wrapper->get_sizes_OZ();
    for (int i = coords[2] - 1; coords[2] > 0; --i) {
      shift_relative_to_other_processes +=
          all_local_Nz[mpi_wrapper->get_cart_rank(coords[0], coords[1], i)];
    }
    int z = bounds.az + steps.dz * shift_relative_to_other_processes;
    int sign = get_sign(E, B);
    int Nz = grid->get_mpi_local_unbounded_subdomain_sizes().Nz;

    for (int i = 0; i < Nx; ++i, z += steps.dz)
      for (int j = 0; j < Ny; ++j)
        for (int k = 0; k < Nz; ++k) {
          int index = Nx * Ny * (k + 1) + Nx * (j + 1) + (i + 1);
          Ex[index] = Ey[index] = Ez[index] = 0.0;
          Bx[index] = By[index] = Bz[index] = 0.0;

          E_field[index] = sin(2.0 * PI * (z - bounds.az - sign * C * t) /
                               (bounds.bz - bounds.az));
          B_field[index] =
              sin(2.0 * PI * (z + steps.dz * coeff - bounds.az - sign * C * t) /
                  (bounds.bz - bounds.az));
        }
  }
  // boundary_synchronization_3D()
} // namespace FDTD
