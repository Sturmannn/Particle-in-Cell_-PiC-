#include "FDTD_MPI.hpp"

namespace FDTD {
FDTD_MPI::FDTD_MPI(std::shared_ptr<Grid> _grid,
                   std::shared_ptr<MPI_Wrapper> _mpi_wrapper)
    : grid(std::move(_grid)), mpi_wrapper(std::move(_mpi_wrapper)) {
  int grid_size = grid->mpi_get_total_bounded_subdomain_sizes();
  Ex = Ey = Ez = Bx = By = Bz = Field(grid_size);
}

void FDTD_MPI::Courant_condition_check(const Shift _shift) const {
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

double FDTD_MPI::get_sign(const Component E, const Component B) const {
  // Ox
  if (E == Component::Ey && B == Component::Bz)
    return 1.0;
  else if (E == Component::Ez && B == Component::By)
    return -1.0;
  // Oy
  else if (E == Component::Ez && B == Component::Bx)
    return 1.0;
  else if (E == Component::Ex && B == Component::Bz)
    return -1.0;
  // Oz
  else if (E == Component::Ex && B == Component::By)
    return 1.0;
  else if (E == Component::Ey && B == Component::Bx)
    return -1.0;
  else {
    std::cout << "\nGet_Sign: Error! Wrong Components! "
              << "E = " << E << " B = " << B << std::endl;
    exit(-1);
  }
}

Field &FDTD_MPI::get_E_field(const Component E) {
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

Field &FDTD_MPI::get_B_field(const Component B) {
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

void FDTD_MPI::boundary_synchronization_3D() {
  // TAG 1. First, each process sends its left face to the left neighbor,
  //    and receives the left face from the right neighbor.
  // TAG 2. Each process sends its right face to the right neighbor,
  //    and receives the right face from the left neighbor.
  // TAG 3. Each process sends its back face to the front face of the back
  // (next) process,
  //    and receives the back face from the front neighbor.
  // TAG 4. Each process sends its front face to the back face of the front
  // (next) process,
  //    and receives the front face from the back neighbor.
  // TAG 5. Each process sends its top face downward to the top neighbor,
  //    and receives the top face from the bottom neighbor.
  // TAG 6. Each process sends its bottom face upward to the bottom neighbor,
  //    and receives the bottom face from the top neighbor.

  int left, right, up, down, front, back;
  MPI_Comm cart_comm = mpi_wrapper->get_cart_comm();

  MPI_Cart_shift(mpi_wrapper->get_cart_comm(), 0, 1, &left, &right);
  MPI_Cart_shift(mpi_wrapper->get_cart_comm(), 1, 1, &back, &front);
  MPI_Cart_shift(mpi_wrapper->get_cart_comm(), 2, 1, &down, &up);

  UnboundedGridSizes grid_sizes =
      grid->get_mpi_local_unbounded_subdomain_sizes();
  const int Nx = grid_sizes.Nx;
  const int Ny = grid_sizes.Ny;
  const int Nz = grid_sizes.Nz;

  const int NyNz = Ny * Nz;
  const int NxNz = (Nx + 2) * Nz;
  const int NxNy = (Nx + 2) * (Ny + 2);
  const int stride_yz = (Ny + 2) * (Nz + 2);
  const int stride_z = Nz + 2;

  int index, storage_index;
  // 1. Send and receive left and right faces
  std::vector<double> left_send(Ny * Nz * 6, 0.0);
  std::vector<double> right_receive(Ny * Nz * 6, 0.0);


  #pragma omp parallel for private(storage_index, index)
  for (int y = 0; y < Ny; ++y)
    #pragma omp simd
    for (int z = 0; z < Nz; ++z) {
      storage_index = y * Nz + z;
      // int index = (Ny + 2) * (Nz + 2) * (0 + 1) + (Nz + 2) * (y + 1) + (z + 1);
      index = stride_yz + stride_z * (y + 1) + (z + 1);
      left_send[storage_index] = Ex[index];
      left_send[storage_index + NyNz * 1] = Ey[index];
      left_send[storage_index + NyNz * 2] = Ez[index];
      left_send[storage_index + NyNz * 3] = Bx[index];
      left_send[storage_index + NyNz * 4] = By[index];
      left_send[storage_index + NyNz * 5] = Bz[index];
    }
  mpi_wrapper->mpi_sendrecv(left_send.data(), Ny * Nz * 6, MPI_DOUBLE, left, 1,
                            right_receive.data(), Ny * Nz * 6, MPI_DOUBLE,
                            right, 1, cart_comm, MPI_STATUS_IGNORE);

  #pragma omp parallel for private(storage_index, index)
  for (int y = 0; y < Ny; ++y)
    #pragma omp simd
    for (int z = 0; z < Nz; ++z) {
      storage_index = y * Nz + z;
      // int index = (Ny + 2) * (Nz + 2) * (Nx + 1) + (Nz + 2) * (y + 1) + (z + 1);
      index = stride_yz * (Nx + 1) + stride_z * (y + 1) + (z + 1);
      Ex[index] = right_receive[storage_index];
      Ey[index] = right_receive[storage_index + NyNz * 1];
      Ez[index] = right_receive[storage_index + NyNz * 2];
      Bx[index] = right_receive[storage_index + NyNz * 3];
      By[index] = right_receive[storage_index + NyNz * 4];
      Bz[index] = right_receive[storage_index + NyNz * 5];
    }

  // 2. Send and receive right and left faces
  std::vector<double> &right_send = left_send;
  std::vector<double> &left_receive = right_receive;

  #pragma omp parallel for private(storage_index, index)
  for (int y = 0; y < Ny; ++y)
    #pragma omp simd
    for (int z = 0; z < Nz; ++z) {
      storage_index = y * Nz + z;
      // int index = (Ny + 2) * (Nz + 2) * ((Nx - 1) + 1) + (Nz + 2) * (y + 1) + (z + 1);
      index = stride_yz * Nx + stride_z * (y + 1) + (z + 1);
      right_send[storage_index] = Ex[index];
      right_send[storage_index + NyNz * 1] = Ey[index];
      right_send[storage_index + NyNz * 2] = Ez[index];
      right_send[storage_index + NyNz * 3] = Bx[index];
      right_send[storage_index + NyNz * 4] = By[index];
      right_send[storage_index + NyNz * 5] = Bz[index];
    }

  mpi_wrapper->mpi_sendrecv(right_send.data(), Ny * Nz * 6, MPI_DOUBLE, right,
                            2, left_receive.data(), Ny * Nz * 6, MPI_DOUBLE,
                            left, 2, cart_comm, MPI_STATUS_IGNORE);

  #pragma omp parallel for private(storage_index, index)
  for (int y = 0; y < Ny; ++y)
    #pragma omp simd
    for (int z = 0; z < Nz; ++z) {
      storage_index = y * Nz + z;
      // int index = (Ny + 2) * (Nz + 2) * (-1 + 1) + (Nz + 2) * (y + 1) + (z + 1);
      index = stride_z * (y + 1) + (z + 1);
      Ex[index] = left_receive[storage_index];
      Ey[index] = left_receive[storage_index + NyNz * 1];
      Ez[index] = left_receive[storage_index + NyNz * 2];
      Bx[index] = left_receive[storage_index + NyNz * 3];
      By[index] = left_receive[storage_index + NyNz * 4];
      Bz[index] = left_receive[storage_index + NyNz * 5];
    }

  // 3. Send and receive back and front faces
  std::vector<double> &back_send = left_send;
  std::vector<double> &front_receive = right_receive;
  back_send.resize((Nx + 2) * Nz * 6);
  front_receive.resize((Nx + 2) * Nz * 6);

  #pragma omp parallel for private(storage_index, index)
  for (int x = -1; x < Nx + 1; ++x)
    #pragma omp simd
    for (int z = 0; z < Nz; ++z) {
      storage_index = (x + 1) * Nz + z;
      // int index = (Ny + 2) * (Nz + 2) * (x + 1) + (Nz + 2) * ((Ny - 1) + 1) + (z + 1);
      index = stride_yz * (x + 1) + stride_z * Ny + (z + 1);
      back_send[storage_index] = Ex[index];
      back_send[storage_index + NxNz * 1] = Ey[index];
      back_send[storage_index + NxNz * 2] = Ez[index];
      back_send[storage_index + NxNz * 3] = Bx[index];
      back_send[storage_index + NxNz * 4] = By[index];
      back_send[storage_index + NxNz * 5] = Bz[index];
    }
  mpi_wrapper->mpi_sendrecv(back_send.data(), (Nx + 2) * Nz * 6, MPI_DOUBLE,
                            back, 3, front_receive.data(), (Nx + 2) * Nz * 6,
                            MPI_DOUBLE, front, 3, cart_comm, MPI_STATUS_IGNORE);

  #pragma omp parallel for private(storage_index, index)
  for (int x = -1; x < Nx + 1; ++x)
    #pragma omp simd
    for (int z = 0; z < Nz; ++z) {
      storage_index = (x + 1) * Nz + z;
      // int index = (Ny + 2) * (Nz + 2) * (x + 1) + (Nz + 2) * (-1 + 1) + (z + 1);
      index = stride_yz * (x + 1) + (z + 1);
      Ex[index] = front_receive[storage_index];
      Ey[index] = front_receive[storage_index + NxNz * 1];
      Ez[index] = front_receive[storage_index + NxNz * 2];
      Bx[index] = front_receive[storage_index + NxNz * 3];
      By[index] = front_receive[storage_index + NxNz * 4];
      Bz[index] = front_receive[storage_index + NxNz * 5];
    }

  // 4. Send and receive front and back faces
  std::vector<double> &front_send = left_send;
  std::vector<double> &back_receive = right_receive;

  #pragma omp parallel for private(storage_index, index)
  for (int x = -1; x < Nx + 1; ++x)
    #pragma omp simd
    for (int z = 0; z < Nz; ++z) {
      storage_index = (x + 1) * Nz + z;
      // int index = (Ny + 2) * (Nz + 2) * (x + 1) + (Nz + 2) * (0 + 1) + (z + 1);
      index = stride_yz * (x + 1) + stride_z + (z + 1);
      front_send[storage_index] = Ex[index];
      front_send[storage_index + NxNz * 1] = Ey[index];
      front_send[storage_index + NxNz * 2] = Ez[index];
      front_send[storage_index + NxNz * 3] = Bx[index];
      front_send[storage_index + NxNz * 4] = By[index];
      front_send[storage_index + NxNz * 5] = Bz[index];
    }

  mpi_wrapper->mpi_sendrecv(front_send.data(), (Nx + 2) * Nz * 6, MPI_DOUBLE,
                            front, 4, back_receive.data(), (Nx + 2) * Nz * 6,
                            MPI_DOUBLE, back, 4, cart_comm, MPI_STATUS_IGNORE);

  #pragma omp parallel for private(storage_index, index)
  for (int x = -1; x < Nx + 1; ++x)
    #pragma omp simd
    for (int z = 0; z < Nz; ++z) {
      storage_index = (x + 1) * Nz + z;
      // int index = (Ny + 2) * (Nz + 2) * (x + 1) + (Nz + 2) * (Ny + 1) + (z + 1);
      index = stride_yz * (x + 1) + stride_z * (Ny + 1) + (z + 1);
      Ex[index] = back_receive[storage_index];
      Ey[index] = back_receive[storage_index + NxNz * 1];
      Ez[index] = back_receive[storage_index + NxNz * 2];
      Bx[index] = back_receive[storage_index + NxNz * 3];
      By[index] = back_receive[storage_index + NxNz * 4];
      Bz[index] = back_receive[storage_index + NxNz * 5];
    }

  // 5. Send and receive top and bottom faces
  std::vector<double> &up_send = left_send;
  std::vector<double> &down_receive = right_receive;
  up_send.resize((Nx + 2) * (Ny + 2) * 6);
  down_receive.resize((Nx + 2) * (Ny + 2) * 6);

  #pragma omp parallel for private(storage_index, index)
  for (int x = -1; x < Nx + 1; ++x)
    #pragma omp simd
    for (int y = -1; y < Ny + 1; ++y) {
      storage_index = (x + 1) * (Ny + 2) + (y + 1);
      // int index = (Ny + 2) * (Nz + 2) * (x + 1) + (Nz + 2) * (y + 1) + ((Nz - 1) + 1);
      index = stride_yz * (x + 1) + stride_z * (y + 1) + Nz;
      up_send[storage_index] = Ex[index];
      up_send[storage_index + NxNy * 1] = Ey[index];
      up_send[storage_index + NxNy * 2] = Ez[index];
      up_send[storage_index + NxNy * 3] = Bx[index];
      up_send[storage_index + NxNy * 4] = By[index];
      up_send[storage_index + NxNy * 5] = Bz[index];
    }

  mpi_wrapper->mpi_sendrecv(up_send.data(), (Nx + 2) * (Ny + 2) * 6, MPI_DOUBLE,
                            up, 5, down_receive.data(), (Nx + 2) * (Ny + 2) * 6,
                            MPI_DOUBLE, down, 5, cart_comm, MPI_STATUS_IGNORE);

  #pragma omp parallel for private(storage_index, index)
  for (int x = -1; x < Nx + 1; ++x)
    #pragma omp simd
    for (int y = -1; y < Ny + 1; ++y) {
      storage_index = (x + 1) * (Ny + 2) + (y + 1);
      // int index = (Ny + 2) * (Nz + 2) * (x + 1) + (Nz + 2) * (y + 1) + (-1 + 1);
      index = stride_yz * (x + 1) + stride_z * (y + 1);
      Ex[index] = down_receive[storage_index];
      Ey[index] = down_receive[storage_index + NxNy * 1];
      Ez[index] = down_receive[storage_index + NxNy * 2];
      Bx[index] = down_receive[storage_index + NxNy * 3];
      By[index] = down_receive[storage_index + NxNy * 4];
      Bz[index] = down_receive[storage_index + NxNy * 5];
    }

  // 6. Send and receive bottom and top faces
  std::vector<double> &down_send = left_send;
  std::vector<double> &up_receive = right_receive;

  #pragma omp parallel for private(storage_index, index)
  for (int x = -1; x < Nx + 1; ++x)
    #pragma omp simd
    for (int y = -1; y < Ny + 1; ++y) {
      storage_index = (x + 1) * (Ny + 2) + (y + 1);
      // int index = (Ny + 2) * (Nz + 2) * (x + 1) + (Nz + 2) * (y + 1) + (0 + 1);
      index = stride_yz * (x + 1) + stride_z * (y + 1) + 1;
      down_send[storage_index] = Ex[index];
      down_send[storage_index + NxNy * 1] = Ey[index];
      down_send[storage_index + NxNy * 2] = Ez[index];
      down_send[storage_index + NxNy * 3] = Bx[index];
      down_send[storage_index + NxNy * 4] = By[index];
      down_send[storage_index + NxNy * 5] = Bz[index];
    }

  mpi_wrapper->mpi_sendrecv(down_send.data(), (Nx + 2) * (Ny + 2) * 6,
                            MPI_DOUBLE, down, 6, up_receive.data(),
                            (Nx + 2) * (Ny + 2) * 6, MPI_DOUBLE, up, 6,
                            cart_comm, MPI_STATUS_IGNORE);

  #pragma omp parallel for private(storage_index, index)
  for (int x = -1; x < Nx + 1; ++x)
    #pragma omp simd
    for (int y = -1; y < Ny + 1; ++y) {
      storage_index = (x + 1) * (Ny + 2) + (y + 1);
      // int index = (Ny + 2) * (Nz + 2) * (x + 1) + (Nz + 2) * (y + 1) + (Nz + 1);
      index = stride_yz * (x + 1) + stride_z * (y + 1) + (Nz + 1);
      Ex[index] = up_receive[storage_index];
      Ey[index] = up_receive[storage_index + NxNy * 1];
      Ez[index] = up_receive[storage_index + NxNy * 2];
      Bx[index] = up_receive[storage_index + NxNy * 3];
      By[index] = up_receive[storage_index + NxNy * 4];
      Bz[index] = up_receive[storage_index + NxNy * 5];
    }
}

Axis FDTD_MPI::get_axis(const Component E, const Component B) {
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

std::string FDTD_MPI::axisToString(const Component E, const Component B) {
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

double FDTD_MPI::get_space_delta(const Component E, const Component B) const {
  Axis axis = get_axis(E, B);
  GridCoordinatesSteps steps = grid->get_coordinates_steps();
  switch (axis) {
  case Axis::X:
    return steps.dx;
  case Axis::Y:
    return steps.dy;
  case Axis::Z:
    return steps.dz;
  default:
    std::cout << "\nGet delta space: Error! Wrong axis!\n";
    exit(-1);
  }
  return -1.0; // Error code
}

void AnalyticalSolverFDTD::solve(const Component E, const Component B,
                                 const double t, const Shift _shift) {
  if (mpi_wrapper->get_world_rank() == 0) {
    std::cout << "Axis: " << axisToString(E, B) << std::endl;
  }

  analytical_soulution(E, B, t, _shift);
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

    for (int i = coords[0]; i > 0; --i) {
      shift_relative_to_other_processes +=
          all_local_Nx[mpi_wrapper->get_cart_rank(i - 1, coords[1], coords[2])];
    }
    double sign = get_sign(E, B);

    // #pragma omp parallel for collapse(3)
    for (int i = 0; i < Nx; ++i) {
      double x = bounds.ax + steps.dx * (shift_relative_to_other_processes + i);
      for (int j = 0; j < Ny; ++j)
        for (int k = 0; k < Nz; ++k) {
          int index =
              (Ny + 2) * (Nz + 2) * (i + 1) + (Nz + 2) * (j + 1) + (k + 1);
          Ex[index] = Ey[index] = Ez[index] = 0.0;
          Bx[index] = By[index] = Bz[index] = 0.0;

          E_field[index] = sin(2.0 * PI * (x - bounds.ax - sign * C * t) /
                               (bounds.bx - bounds.ax));
          B_field[index] =
              sin(2.0 * PI * (x + steps.dx * coeff - bounds.ax - sign * C * t) /
                  (bounds.bx - bounds.ax));
        }
    }
  } else if (axis == Axis::Y) {
    std::vector<int> all_local_Ny = mpi_wrapper->get_sizes_OY();
    for (int j = coords[1]; j > 0; --j) {
      shift_relative_to_other_processes +=
          all_local_Ny[mpi_wrapper->get_cart_rank(coords[0], j - 1, coords[2])];
    }
    double sign = get_sign(E, B);

    // #pragma omp parallel for collapse(3)
    for (int j = 0; j < Ny; ++j) {
      double y = bounds.ay + steps.dy * (shift_relative_to_other_processes + j);
      for (int i = 0; i < Nx; ++i)
        for (int k = 0; k < Nz; ++k) {
          int index =
              (Ny + 2) * (Nz + 2) * (i + 1) + (Nz + 2) * (j + 1) + (k + 1);

          Ex[index] = Ey[index] = Ez[index] = 0.0;
          Bx[index] = By[index] = Bz[index] = 0.0;

          E_field[index] = sin(2.0 * PI * (y - bounds.ay - sign * C * t) /
                               (bounds.by - bounds.ay));

          B_field[index] =
              sin(2.0 * PI * (y + steps.dy * coeff - bounds.ay - sign * C * t) /
                  (bounds.by - bounds.ay));
        }
    }
  } else if (axis == Axis::Z) {
    std::vector<int> all_local_Nz = mpi_wrapper->get_sizes_OZ();
    for (int k = coords[2]; k > 0; --k) {
      shift_relative_to_other_processes +=
          all_local_Nz[mpi_wrapper->get_cart_rank(coords[0], coords[1], k - 1)];
    }
    double sign = get_sign(E, B);

    // #pragma omp parallel for collapse(3)
    for (int k = 0; k < Nz; ++k) {
      double z = bounds.az + steps.dz * (shift_relative_to_other_processes + k);
      for (int j = 0; j < Ny; ++j)
        for (int i = 0; i < Nx; ++i) {
          int index =
              (Ny + 2) * (Nz + 2) * (i + 1) + (Nz + 2) * (j + 1) + (k + 1);
          Ex[index] = Ey[index] = Ez[index] = 0.0;
          Bx[index] = By[index] = Bz[index] = 0.0;

          E_field[index] = sin(2.0 * PI * (z - bounds.az - sign * C * t) /
                               (bounds.bz - bounds.az));
          B_field[index] =
              sin(2.0 * PI * (z + steps.dz * coeff - bounds.az - sign * C * t) /
                  (bounds.bz - bounds.az));
        }
    }
  }
  boundary_synchronization_3D();
}

void NumericalSolverFDTD::solve(const Component E, const Component B,
                                const double t, const Shift _shift) {
  set_default_values(E, B, _shift);
  numerical_solution(t);
}

void NumericalSolverFDTD::update_E_field() {
  const int Nx = grid->get_mpi_local_unbounded_subdomain_sizes().Nx;
  const int Ny = grid->get_mpi_local_unbounded_subdomain_sizes().Ny;
  const int Nz = grid->get_mpi_local_unbounded_subdomain_sizes().Nz;

  // E_dt = grid->get_dt();
  const double C_mul_Edt = C * grid->get_dt();

  const double factor_x = C_mul_Edt / grid->get_coordinates_steps().dx;
  const double factor_y = C_mul_Edt / grid->get_coordinates_steps().dy;
  const double factor_z = C_mul_Edt / grid->get_coordinates_steps().dz;

  const int Ny_mul_Nz = (Ny + 2) * (Nz + 2);

#pragma omp parallel
  {
    int index, index_i1, index_j1, index_k1;
#pragma omp for collapse(2)
    for (int i = 0; i < Nx; ++i)
      for (int j = 0; j < Ny; ++j) {
        index = Ny_mul_Nz * (i + 1) + (Nz + 2) * (j + 1);
#pragma omp simd
        for (int k = 0; k < Nz; ++k) {
          ++index;
          index_i1 = index - Ny_mul_Nz; // i - 1
          index_j1 = index - (Nz + 2);  // j - 1
          index_k1 = index - 1;         // k - 1
          Ex[index] += factor_y * (Bz[index] - Bz[index_j1]) -
                       factor_z * (By[index] - By[index_k1]);
          Ey[index] += factor_z * (Bx[index] - Bx[index_k1]) -
                       factor_x * (Bz[index] - Bz[index_i1]);
          Ez[index] += factor_x * (By[index] - By[index_i1]) -
                       factor_y * (Bx[index] - Bx[index_j1]);
        }
      }
  }
  boundary_synchronization_3D();
}

void NumericalSolverFDTD::update_B_field() {
  const int Nx = grid->get_mpi_local_unbounded_subdomain_sizes().Nx;
  const int Ny = grid->get_mpi_local_unbounded_subdomain_sizes().Ny;
  const int Nz = grid->get_mpi_local_unbounded_subdomain_sizes().Nz;

  // E_dt = grid->get_dt();
  // B_dt = E_dt * 0.5;
  const double C_mul_Bdt = C * grid->get_dt() * 0.5;

  const double factor_x = C_mul_Bdt / grid->get_coordinates_steps().dx;
  const double factor_y = C_mul_Bdt / grid->get_coordinates_steps().dy;
  const double factor_z = C_mul_Bdt / grid->get_coordinates_steps().dz;

  const int Ny_mul_Nz = (Ny + 2) * (Nz + 2);
#pragma omp parallel
  {
    int index, index_i1, index_j1, index_k1;
#pragma omp for collapse(2)
    for (int i = 0; i < Nx; ++i)
      for (int j = 0; j < Ny; ++j) {
        index = Ny_mul_Nz * (i + 1) + (Nz + 2) * (j + 1);
#pragma omp simd
        for (int k = 0; k < Nz; ++k) {
          ++index;
          index_i1 = index + Ny_mul_Nz; // i + 1
          index_j1 = index + (Nz + 2);  // j + 1
          index_k1 = index + 1;         // k + 1

          Bx[index] += factor_z * (Ey[index_k1] - Ey[index]) -
                       factor_y * (Ez[index_j1] - Ez[index]);
          By[index] += factor_x * (Ez[index_i1] - Ez[index]) -
                       factor_z * (Ex[index_k1] - Ex[index]);
          Bz[index] += factor_x * (Ex[index_j1] - Ex[index]) -
                       factor_z * (Ey[index_i1] - Ey[index]);
        }
      }
  }
  boundary_synchronization_3D();
}

void NumericalSolverFDTD::set_default_values(const Component E,
                                             const Component B,
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

    for (int i = coords[0]; i > 0; --i) {
      shift_relative_to_other_processes +=
          all_local_Nx[mpi_wrapper->get_cart_rank(i - 1, coords[1], coords[2])];
    }
    double sign = get_sign(E, B);

    // #pragma omp parallel for collapse(3)
    for (int i = 0; i < Nx; ++i) {
      double x = bounds.ax + steps.dx * (shift_relative_to_other_processes + i);
      for (int j = 0; j < Ny; ++j)
        for (int k = 0; k < Nz; ++k) {
          int index =
              (Ny + 2) * (Nz + 2) * (i + 1) + (Nz + 2) * (j + 1) + (k + 1);
          Ex[index] = Ey[index] = Ez[index] = 0.0;
          Bx[index] = By[index] = Bz[index] = 0.0;

          E_field[index] =
              sin(2.0 * PI * (x - bounds.ax) / (bounds.bx - bounds.ax));
          B_field[index] = sin(2.0 * PI * (x + steps.dx * coeff - bounds.ax) /
                               (bounds.bx - bounds.ax));
        }
    }
  } else if (axis == Axis::Y) {
    std::vector<int> all_local_Ny = mpi_wrapper->get_sizes_OY();
    for (int j = coords[1]; j > 0; --j) {
      shift_relative_to_other_processes +=
          all_local_Ny[mpi_wrapper->get_cart_rank(coords[0], j - 1, coords[2])];
    }
    double sign = get_sign(E, B);

    // #pragma omp parallel for collapse(3)
    for (int j = 0; j < Ny; ++j) {
      double y = bounds.ay + steps.dy * (shift_relative_to_other_processes + j);
      for (int i = 0; i < Nx; ++i)
        for (int k = 0; k < Nz; ++k) {
          int index =
              (Ny + 2) * (Nz + 2) * (i + 1) + (Nz + 2) * (j + 1) + (k + 1);

          Ex[index] = Ey[index] = Ez[index] = 0.0;
          Bx[index] = By[index] = Bz[index] = 0.0;

          E_field[index] =
              sin(2.0 * PI * (y - bounds.ay) / (bounds.by - bounds.ay));

          B_field[index] = sin(2.0 * PI * (y + steps.dy * coeff - bounds.ay) /
                               (bounds.by - bounds.ay));
        }
    }
  } else if (axis == Axis::Z) {
    std::vector<int> all_local_Nz = mpi_wrapper->get_sizes_OZ();
    for (int k = coords[2]; k > 0; --k) {
      shift_relative_to_other_processes +=
          all_local_Nz[mpi_wrapper->get_cart_rank(coords[0], coords[1], k - 1)];
    }
    double sign = get_sign(E, B);

    // #pragma omp parallel for collapse(3)
    for (int k = 0; k < Nz; ++k) {
      double z = bounds.az + steps.dz * (shift_relative_to_other_processes + k);
      for (int j = 0; j < Ny; ++j)
        for (int i = 0; i < Nx; ++i) {
          int index =
              (Ny + 2) * (Nz + 2) * (i + 1) + (Nz + 2) * (j + 1) + (k + 1);
          Ex[index] = Ey[index] = Ez[index] = 0.0;
          Bx[index] = By[index] = Bz[index] = 0.0;

          E_field[index] =
              sin(2.0 * PI * (z - bounds.az) / (bounds.bz - bounds.az));
          B_field[index] = sin(2.0 * PI * (z + steps.dz * coeff - bounds.az) /
                               (bounds.bz - bounds.az));
        }
    }
  }
  boundary_synchronization_3D();
}

void FDTD::NumericalSolverFDTD::numerical_solution(const double t) {
  int iterations = static_cast<int>(t);
  for (int time = 0; time < iterations; ++time) {
    // std::cout << "Time step: " << time << std::endl;
    update_B_field();
    update_E_field();
    update_B_field();
  }
}
} // namespace FDTD
