#include "FDTD.hpp"
#include <string.h>

FDTD::FDTD FDTD::FDTD::create_trash() { 
  return FDTD{{1, 1, 1}, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, 1.0};
}

FDTD::FDTD::FDTD(const std::tuple<int64_t, int64_t, int64_t> &Nx_Ny_Nz,
                 const std::tuple<double, double, double> &ax_ay_az,
                 const std::tuple<double, double, double> &bx_by_bz,
                 const std::tuple<double, double, double> &dx_dy_dz, double _dt)
    : Nx{std::get<0>(Nx_Ny_Nz)}, Ny{std::get<1>(Nx_Ny_Nz)},
      Nz{std::get<2>(Nx_Ny_Nz)}, Ex{get_Nx(), get_Ny(), get_Nz()},
      Ey{get_Nx(), get_Ny(), get_Nz()}, Ez{get_Nx(), get_Ny(), get_Nz()},
      Bx{get_Nx(), get_Ny(), get_Nz()}, By{get_Nx(), get_Ny(), get_Nz()},
      Bz{get_Nx(), get_Ny(), get_Nz()}, ax{std::get<0>(ax_ay_az)},
      ay{std::get<1>(ax_ay_az)}, az{std::get<2>(ax_ay_az)},
      bx{std::get<0>(bx_by_bz)}, by{std::get<1>(bx_by_bz)},
      bz{std::get<2>(bx_by_bz)},
      // dx{(get_bx() - get_ax()) / static_cast<double>(Nx)},
      // dy{(get_by() - get_ay()) / static_cast<double>(Ny)},
      // dz{(Nz > 1) ? (get_bz() - get_az()) / static_cast<double>(Nz) : 0.0},
      dx{std::get<0>(dx_dy_dz)}, dy{std::get<1>(dx_dy_dz)},
      dz{std::get<2>(dx_dy_dz)}, dt{_dt} {
  // int world_size = 0;
  // MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // if (world_size == 1) return;

  // // Определение общего размеры сетки для правильного определения dx, dy, dz

  // int rank = 0;
  // MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // int64_t global_Nx_Ny_Nz[3] = {0, 0, 0};
  // int64_t local_Nx_Ny_Nz[3] = {Nx, Ny, Nz};

  // if (rank == 0) {
  //   std::cout << "Nx = " << local_Nx_Ny_Nz[0] << " Ny = " <<
  //   local_Nx_Ny_Nz[1] << " Nz = " << local_Nx_Ny_Nz[2] << '\n';
  // }

  // MPI_Allreduce(&local_Nx_Ny_Nz, &global_Nx_Ny_Nz, 3, MPI_INT64_T, MPI_SUM,
  // MPI_COMM_WORLD);

  // if (rank == 0) {
  //   std::cout << "Nx = " << global_Nx_Ny_Nz[0] << " Ny = " <<
  //   global_Nx_Ny_Nz[1] << " Nz = " << global_Nx_Ny_Nz[2] << '\n';
  // }

  // global_Nx_Ny_Nz[0] = std::max(global_Nx_Ny_Nz[0] - 2 * world_size,
  // static_cast<int64_t>(1)); global_Nx_Ny_Nz[1] = std::max(global_Nx_Ny_Nz[1]
  // - 2 * world_size, static_cast<int64_t>(1));
  // // Важно правильно настроить Nz, в случае, если сетка 2D
  // // ОПАСНОЕ МЕСТО. Не уверен, учтены ли здесь все случаи!!!
  // global_Nx_Ny_Nz[2] = std::max(global_Nx_Ny_Nz[2] - 2 * world_size,
  // static_cast<int64_t>(1));

  // if (rank == 0) {
  //   std::cout << "Nx = " << global_Nx_Ny_Nz[0] << " Ny = " <<
  //   global_Nx_Ny_Nz[1] << " Nz = " << global_Nx_Ny_Nz[2] << '\n';
  // }

  // dx = (get_bx() - get_ax()) / static_cast<double>(global_Nx_Ny_Nz[0]);
  // dy = (get_by() - get_ay()) / static_cast<double>(global_Nx_Ny_Nz[1]);
  // // Важно правильно настроить dz, в случае, если сетка 2D
  // dz = (global_Nx_Ny_Nz[2] > 1) ? (get_bz() - get_az()) /
  // static_cast<double>(global_Nx_Ny_Nz[2]) : 0.0;

  // std::cout << "rank = " << rank << " dx = " << dx << " dy = " << dy << " dz
  // = " << dz << '\n';
}

FDTD::FDTD::FDTD(const FDTD &_fields)
    : Nx{_fields.get_Nx()}, Ny{_fields.get_Ny()}, Nz{_fields.get_Nz()},

      Ex{_fields.get_Ex()}, Ey{_fields.get_Ey()}, Ez{_fields.get_Ez()},
      Bx{_fields.get_Bx()}, By{_fields.get_By()}, Bz{_fields.get_Bz()},

      ax{_fields.get_ax()}, ay{_fields.get_ay()}, az{_fields.get_az()},

      bx{_fields.get_bx()}, by{_fields.get_by()}, bz{_fields.get_bz()},

      dx{_fields.get_dx()}, dy{_fields.get_dy()}, dz{_fields.get_dz()},

      dt{_fields.get_dt()} {}

FDTD::FDTD::FDTD(FDTD &&_fields) noexcept
    : Nx{_fields.get_Nx()}, Ny{_fields.get_Ny()}, Nz{_fields.get_Nz()},
      Ex{std::move(_fields.get_Ex())}, Ey{std::move(_fields.get_Ey())},
      Ez{std::move(_fields.get_Ez())}, Bx{std::move(_fields.get_Bx())},
      By{std::move(_fields.get_By())}, Bz{std::move(_fields.get_Bz())},
      ax{_fields.get_ax()}, ay{_fields.get_ay()}, az{_fields.get_az()},
      bx{_fields.get_bx()}, by{_fields.get_by()}, bz{_fields.get_bz()},
      dx{_fields.get_dx()}, dy{_fields.get_dy()}, dz{_fields.get_dz()},
      dt{_fields.get_dt()} {}

FDTD::FDTD &FDTD::FDTD::operator=(const FDTD &_fields) {
  if (this != &_fields) {
    // this->FDTD::FDTD(_fields);
    Nx = _fields.get_Nx();
    Ny = _fields.get_Ny();
    Nz = _fields.get_Nz();
    Ex = _fields.get_Ex();
    Ey = _fields.get_Ey();
    Ez = _fields.get_Ez();
    Bx = _fields.get_Bx();
    By = _fields.get_By();
    Bz = _fields.get_Bz();
    ax = _fields.get_ax();
    ay = _fields.get_ay();
    az = _fields.get_az();
    bx = _fields.get_bx();
    by = _fields.get_by();
    bz = _fields.get_bz();
    dx = _fields.get_dx();
    dy = _fields.get_dy();
    dz = _fields.get_dz();
    dt = _fields.get_dt();
  }
  return *this;
}

FDTD::FDTD &FDTD::FDTD::operator=(FDTD &&_fields) noexcept {
  if (this != &_fields) {
    // this->FDTD::FDTD(std::move(_fields));
    Nx = _fields.get_Nx();
    Ny = _fields.get_Ny();
    Nz = _fields.get_Nz();
    Ex = std::move(_fields.get_Ex());
    Ey = std::move(_fields.get_Ey());
    Ez = std::move(_fields.get_Ez());
    Bx = std::move(_fields.get_Bx());
    By = std::move(_fields.get_By());
    Bz = std::move(_fields.get_Bz());
    ax = _fields.get_ax();
    ay = _fields.get_ay();
    az = _fields.get_az();
    bx = _fields.get_bx();
    by = _fields.get_by();
    bz = _fields.get_bz();
    dx = _fields.get_dx();
    dy = _fields.get_dy();
    dz = _fields.get_dz();
    dt = _fields.get_dt();
  }
  return *this;
}

void FDTD::FDTD::set_Nx_Ny_Nz(int64_t _Nx, int64_t _Ny, int64_t _Nz) {
  Nx = _Nx;
  Ny = _Ny;

  Nz == 0 ? Nz = 1 : Nz = _Nz;

  dx = (bx - ax) / static_cast<double>(Nx);
  dy = (by - ay) / static_cast<double>(Ny);
  dz = (bz - az) / static_cast<double>(Nz);

  Ex.resize_field(_Nx, _Ny, _Nz);
  Ey.resize_field(_Nx, _Ny, _Nz);
  Ez.resize_field(_Nx, _Ny, _Nz);
  Bx.resize_field(_Nx, _Ny, _Nz);
  By.resize_field(_Nx, _Ny, _Nz);
  Bz.resize_field(_Nx, _Ny, _Nz);
}

// Only for 2D
void FDTD::FDTD::field_update(const double t, MPI_Comm cart_comm) {
  if (dt == 0.0) {
    std::cout << "Time step is null";
    exit(-1);
  }

  std::cout << "unshifted_field_update(const double t)" << std::endl;

  int64_t i{0};
  int64_t j{0};

#pragma omp parallel private(i, j)
{
  for (double time = 0.0; time < t; time += dt) {
#pragma omp for collapse(2)
    for (j = 0; j < Ny; ++j) {
      for (i = 0; i < Nx; ++i) {
        Ex(i, j) =
            Ex(i, j) + C * dt * 0.5 * ((Bz(i, j + 1) - Bz(i, j - 1)) / dy);
        Ey(i, j) =
            Ey(i, j) - C * dt * 0.5 * ((Bz(i + 1, j) - Bz(i - 1, j)) / dx);
        Ez(i, j) = Ez(i, j) + C * dt * 0.5 *
                                  (((By(i + 1, j) - By(i - 1, j)) / dx) -
                                   (Bx(i, j + 1) - Bx(i, j - 1)) / dy);
      }
    }
    boundary_synchronization(cart_comm);
#pragma omp for collapse(2)
    for (j = 0; j < Ny; ++j)
      for (i = 0; i < Nx; ++i) {
        Bx(i, j) =
            Bx(i, j) - C * dt * 0.5 * ((Ez(i, j + 1) - Ez(i, j - 1)) / dy);
        By(i, j) =
            By(i, j) + C * dt * 0.5 * ((Ez(i + 1, j) - Ez(i - 1, j)) / dx);
        Bz(i, j) = Bz(i, j) - C * dt * 0.5 *
                                  (((Ey(i + 1, j) - Ey(i - 1, j)) / dx) -
                                   (Ex(i, j + 1) - Ex(i, j - 1)) / dy);
      }
    boundary_synchronization(cart_comm);
  }
}
}

// Only for 2D
void FDTD::FDTD::field_update(const int64_t t, MPI_Comm cart_comm) {
  if (dt == 0.0) {
    std::cout << "Time step is null";
    exit(-1);
  }

  std::cout << "unshifted_field_update(const int64_t t)" << std::endl;

  int64_t i{0};
  int64_t j{0};

#pragma omp parallel private(i, j)
{
  for (int64_t time = 0; time < t; time++) {
#pragma omp for collapse(2)
    for (j = 0; j < Ny; ++j) {
      for (i = 0; i < Nx; ++i) {
        Ex(i, j) =
            Ex(i, j) + C * dt * 0.5 * ((Bz(i, j + 1) - Bz(i, j - 1)) / dy);
        Ey(i, j) =
            Ey(i, j) - C * dt * 0.5 * ((Bz(i + 1, j) - Bz(i - 1, j)) / dx);
        Ez(i, j) = Ez(i, j) + C * dt * 0.5 *
                                  (((By(i + 1, j) - By(i - 1, j)) / dx) -
                                   (Bx(i, j + 1) - Bx(i, j - 1)) / dy);
      }
    }
    boundary_synchronization(cart_comm);
#pragma omp for collapse(2)
    for (j = 0; j < Ny; ++j)
      for (i = 0; i < Nx; ++i) {
        Bx(i, j) =
            Bx(i, j) - C * dt * 0.5 * ((Ez(i, j + 1) - Ez(i, j - 1)) / dy);
        By(i, j) =
            By(i, j) + C * dt * 0.5 * ((Ez(i + 1, j) - Ez(i - 1, j)) / dx);
        Bz(i, j) = Bz(i, j) - C * dt * 0.5 *
                                  (((Ey(i + 1, j) - Ey(i - 1, j)) / dx) -
                                   (Ex(i, j + 1) - Ex(i, j - 1)) / dy);
      }
    boundary_synchronization(cart_comm);
  }
}
}

void FDTD::FDTD::shifted_field_update(const double t, MPI_Comm cart_comm) {}
void FDTD::FDTD::shifted_field_update(const int64_t t, MPI_Comm cart_comm) {
  if (dt == 0.0) {
    std::cout << "Time step is null";
    exit(-1);
  }

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
  std::cout << "Start shifted_field_update(const int64_t t)" << std::endl;

  double B_dt = dt * 0.5;
  double E_dt = dt;

  int64_t _Nz = (Nz > 1) ? Nz : 0; // Учитываем случай 2D
  int64_t k_value = (_Nz > 1) ? 0 : -1; // Это эквивалентно -1 для 2D

  int64_t i{0};
  int64_t j{0};
  int64_t k{k_value};

  if (rank == 0) {
    std::cout << "Nz = " << Nz << " _Nz = " << _Nz << " k_value = " << k_value << std::endl;
    std::cout << "dx = " << dx << " dy = " << dy << " dz = " << dz << std::endl;
  }

  // for (double time = 0.0; time < t; time += E_dt) // ПОПРАВИТЬ НА time += E_dt
  if (Nz < 2) {
// #pragma omp parallel private(i, j) 
  {
    for (int64_t time = 0; time < t; ++time) {

// #pragma omp for collapse(2)
      for (i = 0; i < Nx; ++i)
        for (j = 0; j < Ny; ++j) {
          Bx(i, j) = Bx(i, j) + C * B_dt * (-(Ez(i, j + 1) - Ez(i, j)) / dy);
          By(i, j) = By(i, j) + C * B_dt * ((Ez(i + 1, j) - Ez(i, j)) / dx);
          Bz(i, j) = Bz(i, j) + C * B_dt *
                                    ((Ex(i, j + 1) - Ex(i, j)) / dy -
                                     (Ey(i + 1, j) - Ey(i, j)) / dx);

        }
      // for (i = 0; i < Nx; ++i)
      //   for (j = 0; j < Ny; ++j) {
      //     if (Bz(i, j) - Bz(i, j - 1) != 0.0) {
      //       std::cout << "Ex = " << Ex(i,j) << " " << "i = " << i << " j = " << j << " Bz(i , j) - Bz(i, j - 1) = " << Bz(i, j) - Bz(i, j - 1) << " rank = " << rank << " time = " << time << std::endl;
      //       exit(-10);
      //     }      
      //   }                               
      boundary_synchronization(cart_comm);
// #pragma omp for collapse(2)
      for (i = 0; i < Nx; ++i)
        for (j = 0; j < Ny; ++j) {
          Ex(i, j) = Ex(i, j) + C * E_dt * ((Bz(i, j) - Bz(i, j - 1)) / dy);
          Ey(i, j) = Ey(i, j) + C * E_dt * (-(Bz(i, j) - Bz(i - 1, j)) / dx);
          Ez(i, j) = Ez(i, j) + C * E_dt *
                                    ((By(i, j) - By(i - 1, j)) / dx -
                                     (Bx(i, j) - Bx(i, j - 1)) / dy);

        }
      boundary_synchronization(cart_comm);
// #pragma omp for collapse(2)
      for (i = 0; i < Nx; ++i)
        for (j = 0; j < Ny; ++j) {
          // if (rank == 0)
            // std::cout << " time = " << time << " i = " << i << " j = " << j << std::endl;
          Bx(i, j) = Bx(i, j) + C * B_dt * (-(Ez(i, j + 1) - Ez(i, j)) / dy);
          By(i, j) = By(i, j) + C * B_dt * ((Ez(i + 1, j) - Ez(i, j)) / dx);
          Bz(i, j) = Bz(i, j) + C * B_dt *
                                    ((Ex(i, j + 1) - Ex(i, j)) / dy -
                                     (Ey(i + 1, j) - Ey(i, j)) / dx);
        }
      boundary_synchronization(cart_comm);
    }
  }
  }  else {
// #pragma omp parallel private(i, j, k)
    {
      for (int64_t time = 0; time < t; ++time) // ПОПРАВИТЬ НА time += E_dt
      {
        // std::cout << "Hello from " << omp_get_thread_num() << std::endl;
// #pragma omp for collapse(3)
        for (i = 0; i < Nx; ++i)
          for (j = 0; j < Ny; ++j)
            for (k = k_value; k < _Nz; ++k) {
              Bx(i, j, k) =
                  Bx(i, j, k) + C * B_dt *
                                    ((Ey(i, j, k + 1) - Ey(i, j, k)) / dz -
                                     (Ez(i, j + 1, k) - Ez(i, j, k)) / dy);

              By(i, j, k) =
                  By(i, j, k) + C * B_dt *
                                    ((Ez(i + 1, j, k) - Ez(i, j, k)) / dx -
                                     (Ex(i, j, k + 1) - Ex(i, j, k)) / dz);

              Bz(i, j, k) =
                  Bz(i, j, k) + C * B_dt *
                                    ((Ex(i, j + 1, k) - Ex(i, j, k)) / dy -
                                     (Ey(i + 1, j, k) - Ey(i, j, k)) / dx);
            }
        boundary_synchronization_3D(cart_comm);
// #pragma omp for collapse(3)
        for (i = 0; i < Nx; ++i)
          for (j = 0; j < Ny; ++j)
            for (k = k_value; k < _Nz; ++k) {
              Ex(i, j, k) =
                  Ex(i, j, k) + C * E_dt *
                                    ((Bz(i, j, k) - Bz(i, j - 1, k)) / dy -
                                     (By(i, j, k) - By(i, j, k - 1)) / dz);
              Ey(i, j, k) =
                  Ey(i, j, k) + C * E_dt *
                                    ((Bx(i, j, k) - Bx(i, j, k - 1)) / dz -
                                     (Bz(i, j, k) - Bz(i - 1, j, k)) / dx);
              Ez(i, j, k) =
                  Ez(i, j, k) + C * E_dt *
                                    ((By(i, j, k) - By(i - 1, j, k)) / dx -
                                     (Bx(i, j, k) - Bx(i, j - 1, k)) / dy);
            }
        boundary_synchronization_3D(cart_comm);
// #pragma omp for collapse(3)
        for (i = 0; i < Nx; ++i)
          for (j = 0; j < Ny; ++j)
            for (k = k_value; k < _Nz; ++k) {
              Bx(i, j, k) =
                  Bx(i, j, k) + C * B_dt *
                                    ((Ey(i, j, k + 1) - Ey(i, j, k)) / dz -
                                     (Ez(i, j + 1, k) - Ez(i, j, k)) / dy);
              By(i, j, k) =
                  By(i, j, k) + C * B_dt *
                                    ((Ez(i + 1, j, k) - Ez(i, j, k)) / dx -
                                     (Ex(i, j, k + 1) - Ex(i, j, k)) / dz);
              Bz(i, j, k) =
                  Bz(i, j, k) + C * B_dt *
                                    ((Ex(i, j + 1, k) - Ex(i, j, k)) / dy -
                                     (Ey(i + 1, j, k) - Ey(i, j, k)) / dx);
            }
        boundary_synchronization_3D(cart_comm);
      }
    }
  }

 if (rank == 0) {
    std::cout << "End shifted_field_update(const int64_t t)" << std::endl;
  }
}


void FDTD::FDTD::write_fields_to_file(const char *directory_path, const Component E,
                                      const Component B, const double delta,
                                      const int64_t row_number) {

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  fs::path file_path = fs::path(directory_path);

  // Формирование имени файла для текущего процесса
  // Обращение сначала к patent_path из-за того, что путь заканчивается на / (поэтому последний элемент пустой)
  if (file_path.parent_path().filename() == "my_data")
    file_path /= "my_data_" + std::to_string(rank) + ".csv";
  else if (file_path.parent_path().filename() == "analytical_data")
    file_path /= "analytical_data_" + std::to_string(rank) + ".csv";
  else {
    std::cerr << "FDTD::FDTD::write_fields_to_file() error: invalid directory path" << std::endl;
    exit(-1);
  }

  Axis axis = get_axis(E, B);
  switch (axis) {
  case Axis::Ox:
    Ex.write_field_to_file_OX(file_path.string().c_str(), row_number);
    Ey.write_field_to_file_OX(file_path.string().c_str(), row_number);
    Ez.write_field_to_file_OX(file_path.string().c_str(), row_number);
    Bx.write_field_to_file_OX(file_path.string().c_str(), row_number);
    By.write_field_to_file_OX(file_path.string().c_str(), row_number);
    Bz.write_field_to_file_OX(file_path.string().c_str(), row_number);
    break;
  case Axis::Oy:
    Ex.write_field_to_file_OY(file_path.string().c_str(), row_number);
    Ey.write_field_to_file_OY(file_path.string().c_str(), row_number);
    Ez.write_field_to_file_OY(file_path.string().c_str(), row_number);
    Bx.write_field_to_file_OY(file_path.string().c_str(), row_number);
    By.write_field_to_file_OY(file_path.string().c_str(), row_number);
    Bz.write_field_to_file_OY(file_path.string().c_str(), row_number);
    break;
  case Axis::Oz:
    Ex.write_field_to_file_OZ(file_path.string().c_str(), row_number);
    Ey.write_field_to_file_OZ(file_path.string().c_str(), row_number);
    Ez.write_field_to_file_OZ(file_path.string().c_str(), row_number);
    Bx.write_field_to_file_OZ(file_path.string().c_str(), row_number);
    By.write_field_to_file_OZ(file_path.string().c_str(), row_number);
    Bz.write_field_to_file_OZ(file_path.string().c_str(), row_number);
    break;
  default:
    break;
  }
  std::ofstream outfile;
  outfile.open(file_path, std::ios::app);
  if (!outfile.is_open()) {
    std::cout << "The file can't be opened!" << std::endl;
    exit(-1);
  }
  outfile << delta << std::endl;
  outfile.close();
}

void FDTD::FDTD::clear_fields(void) noexcept {
  Ex.clear_field();
  Ey.clear_field();
  Ez.clear_field();
  Bx.clear_field();
  By.clear_field();
  Bz.clear_field();
  Nx = Ny = Nz = 0;
  ax = ay = az = 0.0;
  bx = by = bz = 0.0;
  dx = dy = dz = 0.0;
  dt = 0.0;
}

Axis FDTD::FDTD::get_axis(const Component E, const Component B) {
  if ((E == Component::Ey && B == Component::Bz) ||
      (E == Component::Ez && B == Component::By))
    return Axis::Ox;
  else if ((E == Component::Ez && B == Component::Bx) ||
           (E == Component::Ex && B == Component::Bz))
    return Axis::Oy;
  else if ((E == Component::Ex && B == Component::By) ||
           (E == Component::Ey && B == Component::Bx))
    return Axis::Oz;
  else {
    std::cout << "\nGet_Axis: Error! Wrong Components! " << "E = " << E << " B = " << B << std::endl;
    exit(-1);
  }
}

std::string FDTD::FDTD::axisToString(const Component E, const Component B) {
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
    std::cout << "\nGet_Axis: Error! Wrong Components! " << "E = " << E << " B = " << B << std::endl;
    exit(-1);
  }
}

// void FDTD::FDTD::boundary_synchronization(/*Сюда добавить поля*/) {
void FDTD::FDTD::boundary_synchronization(MPI_Comm cart_comm) {

  // СИСТЕМА ТЕГОВ ДЛЯ MPI (ОТПРАВКА):
  // 0 - низ
  // 1 - верх

  int rank = 0;
  int mpi_comm_size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_comm_size);


  // Инструкция по синхронизации границ для 2D в случае разбиения по двум осям:
  // Границы отправляются по отдельности
  // ТЕГ 1. Сначала отпраляется верхняя часть себя - в нижнюю часть верхнего соседа, 
  //    и приём верхней части нижнего соседа - в нижнюю часть себя

  // ТЕГ 2. Потом отправляется правая часть себя - в левую часть правого соседа,
  //    и приём правой части левого соседа - влево себе

  // ТЕГ 3. Потом отправляется нижняя часть себя - в верхнюю часть нижнего соседу,
  //    и приём нижней части верхнего соседа - в верхнюю часть себя

  // ТЕГ 4. Потом отправляется левая часть себя - в правую часть левого соседа,
  //    и приём левой части правого соседа - вправо себе

  // После этого обновляются угловые значения (работают диагональные соседи)
  // ТЕГ 5. Верхний левый угол себя - в правый нижний угол левого верхнего соседа,
  //    и приём верхнего левого угла правого нижнего соседа - в правый нижний угол себя

  // ТЕГ 6. Верхний правый угол себя - в левый нижний угол правого верхнего соседа,
  //    и приём верхнего правого угла левого нижнего соседа - в левый нижний угол себя

  // ТЕГ 7. Нижний правый угол себя - в левый верхний угол правого нижнего соседа,
  //    и приём нижнего правого угла левого верхнего соседа - в левый верхний угол себя

  // ТЕГ 8. Нижний левый угол себя - в правый верхний угол левого нижнего соседа,
  //    и приём нижнего левого угла правого верхнего соседа - в правый верхний угол себя

// ВОПРОС: КАК БУДЕТ РАБОТАТЬ ДЕКАРТОВА ТОПОЛОГИЯ ПРИ 1 ПРОЦЕССЕ?!?!??!?!
// Также пока что рассматривается 2д версия, поэтому третий индекс -1, который обращается в 0. (ДОЛЖЕН БЫТЬ -1??)

// НА ДАННЫЙ МОМЕНТ ПОКА НЕПОНЯТНО, ЧЁ С ИНДЕКСАЦИЕЙ. СДВИГИ ТО ПОМЕНЯНЫ, НО ВНУТРИ ПЕРЕДАЧ С coords НАДО МЕНЯТЬ

// Получение рангов соседних процессов
int left, right, up, down;
MPI_Cart_shift(cart_comm, 0, 1, &left, &right);
MPI_Cart_shift(cart_comm, 1, 1, &up, &down);

// MPI_Cart_shift(cart_comm, 0, 1, &up, &down);
// MPI_Cart_shift(cart_comm, 1, 1, &left, &right);

// MPI_Barrier(cart_comm);
// std::cout << "rank = " << rank << " left = " << left << " right = " << right << " up = " << up << " down = " << down << std::endl;

if (left == MPI_PROC_NULL || right == MPI_PROC_NULL || up == MPI_PROC_NULL || down == MPI_PROC_NULL) {
  std::cout << "FDTD::FDTD::boundary_synchronization() neighbor rank error: MPI_PROC_NULL" << std::endl;
  exit(-1);
}



// 2 ПУНКТ ---------------------------------------------------

// Векторы для отправки и приёма боковых элементов
std::vector<std::vector<double>> right_send(6, std::vector<double>(Ny, 0.0));
std::vector<std::vector<double>> left_receive(6, std::vector<double>(Ny, 0.0));

// Заполняем векторы данными для отправки
for (int64_t j = 0; j < Ny; ++j) {
    right_send[0][j] = Ex(Nx - 1, j); // Ex
    right_send[1][j] = Ey(Nx - 1, j); // Ey
    right_send[2][j] = Ez(Nx - 1, j); // Ez
    right_send[3][j] = Bx(Nx - 1, j); // Bx
    right_send[4][j] = By(Nx - 1, j); // By
    right_send[5][j] = Bz(Nx - 1, j); // Bz
}

MPI_Sendrecv(right_send[0].data(), Ny, MPI_DOUBLE, right, 2, left_receive[0].data(), Ny, MPI_DOUBLE, left, 2, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(right_send[1].data(), Ny, MPI_DOUBLE, right, 2, left_receive[1].data(), Ny, MPI_DOUBLE, left, 2, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(right_send[2].data(), Ny, MPI_DOUBLE, right, 2, left_receive[2].data(), Ny, MPI_DOUBLE, left, 2, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(right_send[3].data(), Ny, MPI_DOUBLE, right, 2, left_receive[3].data(), Ny, MPI_DOUBLE, left, 2, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(right_send[4].data(), Ny, MPI_DOUBLE, right, 2, left_receive[4].data(), Ny, MPI_DOUBLE, left, 2, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(right_send[5].data(), Ny, MPI_DOUBLE, right, 2, left_receive[5].data(), Ny, MPI_DOUBLE, left, 2, cart_comm, MPI_STATUS_IGNORE);

// Записываем данные из вектора приёма в поля
for (int64_t j = 0; j < Ny; ++j) {
    Ex(-1, j) = left_receive[0][j]; // Ex
    Ey(-1, j) = left_receive[1][j]; // Ey
    Ez(-1, j) = left_receive[2][j]; // Ez
    Bx(-1, j) = left_receive[3][j]; // Bx
    By(-1, j) = left_receive[4][j]; // By
    Bz(-1, j) = left_receive[5][j]; // Bz
}

// 4 ПУНКТ ---------------------------------------------------

// Векторы для отправки и приёма боковых элементов
std::vector<std::vector<double>> left_send(6, std::vector<double>(Ny, 0.0));
std::vector<std::vector<double>> right_receive(6, std::vector<double>(Ny, 0.0));

// Заполняем векторы данными для отправки
for (int64_t j = 0; j < Ny; ++j) {
    left_send[0][j] = Ex(0, j); // Ex
    left_send[1][j] = Ey(0, j); // Ey
    left_send[2][j] = Ez(0, j); // Ez
    left_send[3][j] = Bx(0, j); // Bx
    left_send[4][j] = By(0, j); // By
    left_send[5][j] = Bz(0, j); // Bz
}

MPI_Sendrecv(left_send[0].data(), Ny, MPI_DOUBLE, left, 4, right_receive[0].data(), Ny, MPI_DOUBLE, right, 4, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(left_send[1].data(), Ny, MPI_DOUBLE, left, 4, right_receive[1].data(), Ny, MPI_DOUBLE, right, 4, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(left_send[2].data(), Ny, MPI_DOUBLE, left, 4, right_receive[2].data(), Ny, MPI_DOUBLE, right, 4, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(left_send[3].data(), Ny, MPI_DOUBLE, left, 4, right_receive[3].data(), Ny, MPI_DOUBLE, right, 4, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(left_send[4].data(), Ny, MPI_DOUBLE, left, 4, right_receive[4].data(), Ny, MPI_DOUBLE, right, 4, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(left_send[5].data(), Ny, MPI_DOUBLE, left, 4, right_receive[5].data(), Ny, MPI_DOUBLE, right, 4, cart_comm, MPI_STATUS_IGNORE);

// Записываем данные из вектора приёма в поля
for (int64_t j = 0; j < Ny; ++j) {
    Ex(Nx, j) = right_receive[0][j]; // Ex
    Ey(Nx, j) = right_receive[1][j]; // Ey
    Ez(Nx, j) = right_receive[2][j]; // Ez
    Bx(Nx, j) = right_receive[3][j]; // Bx
    By(Nx, j) = right_receive[4][j]; // By
    Bz(Nx, j) = right_receive[5][j]; // Bz
}

// 1 ПУНКТ ---------------------------------------------------

// MPI_Sendrecv(&Ex(0, Ny - 1), Nx, MPI_DOUBLE, up, 1, &Ex(0, -1), Nx, MPI_DOUBLE, down, 1, cart_comm, MPI_STATUS_IGNORE);
// MPI_Sendrecv(&Ey(0, Ny - 1), Nx, MPI_DOUBLE, up, 1, &Ey(0, -1), Nx, MPI_DOUBLE, down, 1, cart_comm, MPI_STATUS_IGNORE);
// MPI_Sendrecv(&Ez(0, Ny - 1), Nx, MPI_DOUBLE, up, 1, &Ez(0, -1), Nx, MPI_DOUBLE, down, 1, cart_comm, MPI_STATUS_IGNORE);
// MPI_Sendrecv(&Bx(0, Ny - 1), Nx, MPI_DOUBLE, up, 1, &Bx(0, -1), Nx, MPI_DOUBLE, down, 1, cart_comm, MPI_STATUS_IGNORE);
// MPI_Sendrecv(&By(0, Ny - 1), Nx, MPI_DOUBLE, up, 1, &By(0, -1), Nx, MPI_DOUBLE, down, 1, cart_comm, MPI_STATUS_IGNORE);
// MPI_Sendrecv(&Bz(0, Ny - 1), Nx, MPI_DOUBLE, up, 1, &Bz(0, -1), Nx, MPI_DOUBLE, down, 1, cart_comm, MPI_STATUS_IGNORE);

MPI_Sendrecv(&Ex(-1, Ny - 1), Nx + 2, MPI_DOUBLE, up, 1, &Ex(-1, -1), Nx + 2, MPI_DOUBLE, down, 1, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(&Ey(-1, Ny - 1), Nx + 2, MPI_DOUBLE, up, 1, &Ey(-1, -1), Nx + 2, MPI_DOUBLE, down, 1, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(&Ez(-1, Ny - 1), Nx + 2, MPI_DOUBLE, up, 1, &Ez(-1, -1), Nx + 2, MPI_DOUBLE, down, 1, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(&Bx(-1, Ny - 1), Nx + 2, MPI_DOUBLE, up, 1, &Bx(-1, -1), Nx + 2, MPI_DOUBLE, down, 1, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(&By(-1, Ny - 1), Nx + 2, MPI_DOUBLE, up, 1, &By(-1, -1), Nx + 2, MPI_DOUBLE, down, 1, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(&Bz(-1, Ny - 1), Nx + 2, MPI_DOUBLE, up, 1, &Bz(-1, -1), Nx + 2, MPI_DOUBLE, down, 1, cart_comm, MPI_STATUS_IGNORE);

// 3 ПУНКТ ---------------------------------------------------

// MPI_Sendrecv(&Ex(0, 0), Nx, MPI_DOUBLE, down, 3, &Ex(0, Ny), Nx, MPI_DOUBLE, up, 3, cart_comm, MPI_STATUS_IGNORE);
// MPI_Sendrecv(&Ey(0, 0), Nx, MPI_DOUBLE, down, 3, &Ey(0, Ny), Nx, MPI_DOUBLE, up, 3, cart_comm, MPI_STATUS_IGNORE);
// MPI_Sendrecv(&Ez(0, 0), Nx, MPI_DOUBLE, down, 3, &Ez(0, Ny), Nx, MPI_DOUBLE, up, 3, cart_comm, MPI_STATUS_IGNORE);
// MPI_Sendrecv(&Bx(0, 0), Nx, MPI_DOUBLE, down, 3, &Bx(0, Ny), Nx, MPI_DOUBLE, up, 3, cart_comm, MPI_STATUS_IGNORE);
// MPI_Sendrecv(&By(0, 0), Nx, MPI_DOUBLE, down, 3, &By(0, Ny), Nx, MPI_DOUBLE, up, 3, cart_comm, MPI_STATUS_IGNORE);
// MPI_Sendrecv(&Bz(0, 0), Nx, MPI_DOUBLE, down, 3, &Bz(0, Ny), Nx, MPI_DOUBLE, up, 3, cart_comm, MPI_STATUS_IGNORE);

MPI_Sendrecv(&Ex(-1, 0), Nx + 2, MPI_DOUBLE, down, 3, &Ex(-1, Ny), Nx + 2, MPI_DOUBLE, up, 3, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(&Ey(-1, 0), Nx + 2, MPI_DOUBLE, down, 3, &Ey(-1, Ny), Nx + 2, MPI_DOUBLE, up, 3, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(&Ez(-1, 0), Nx + 2, MPI_DOUBLE, down, 3, &Ez(-1, Ny), Nx + 2, MPI_DOUBLE, up, 3, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(&Bx(-1, 0), Nx + 2, MPI_DOUBLE, down, 3, &Bx(-1, Ny), Nx + 2, MPI_DOUBLE, up, 3, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(&By(-1, 0), Nx + 2, MPI_DOUBLE, down, 3, &By(-1, Ny), Nx + 2, MPI_DOUBLE, up, 3, cart_comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(&Bz(-1, 0), Nx + 2, MPI_DOUBLE, down, 3, &Bz(-1, Ny), Nx + 2, MPI_DOUBLE, up, 3, cart_comm, MPI_STATUS_IGNORE);



// // 5 ПУНКТ ---------------------------------------------------

// // Получаем количество измерений в топологии
// int ndims;
// MPI_Cartdim_get(cart_comm, &ndims);

// // Получаем размеры топологии, периодичность и координаты текущего процесса
// int dims[ndims];
// int periods[ndims];
// int coords[ndims];
// MPI_Cart_get(cart_comm, ndims, dims, periods, coords);

// // Вычисляем координаты левого верхнего соседа
// // Здесь сразу же используется периодичность (без неё работать не будет)
// int left_upper_coords[2] = {(coords[0] - 1 + dims[0]) % dims[0], (coords[1] + 1) % dims[1]};

// // Получаем ранг левого верхнего соседа
// int left_upper_rank;
// MPI_Cart_rank(cart_comm, left_upper_coords, &left_upper_rank);

// // Вычисляем координаты правого нижнего соседа
// int right_lower_coords[2] = {(coords[0] + 1) % dims[0], (coords[1] - 1 + dims[1]) % dims[1]};

// // Получаем ранг правого нижнего соседа
// int right_lower_rank;
// MPI_Cart_rank(cart_comm, right_lower_coords, &right_lower_rank);

// // Векторы для отправки и приёма угловых элементов
// std::vector<double> left_upper_send(6, 0.0);
// std::vector<double> right_lower_receive(6, 0.0);

// // Заполняем векторы данными для отправки
// left_upper_send[0] = Ex(0, Ny - 1); // Ex
// left_upper_send[1] = Ey(0, Ny - 1); // Ey
// left_upper_send[2] = Ez(0, Ny - 1); // Ez
// left_upper_send[3] = Bx(0, Ny - 1); // Bx
// left_upper_send[4] = By(0, Ny - 1); // By
// left_upper_send[5] = Bz(0, Ny - 1); // Bz

// // Отправляем и получаем данные
// MPI_Sendrecv(left_upper_send.data(), 6, MPI_DOUBLE, left_upper_rank, 5,
//  right_lower_receive.data(), 6, MPI_DOUBLE, right_lower_rank, 5, cart_comm, MPI_STATUS_IGNORE);

// // Записываем данные из вектора приёма в поля
// Ex(Nx, -1) = right_lower_receive[0]; // Ex
// Ey(Nx, -1) = right_lower_receive[1]; // Ey
// Ez(Nx, -1) = right_lower_receive[2]; // Ez
// Bx(Nx, -1) = right_lower_receive[3]; // Bx
// By(Nx, -1) = right_lower_receive[4]; // By
// Bz(Nx, -1) = right_lower_receive[5]; // Bz

// // ПУНКТ 6 ---------------------------------------------------

// // Векторы для отправки и приёма угловых элементов
// std::vector<double> right_upper_send(6, 0.0);
// std::vector<double> left_lower_receive(6, 0.0);

// // Заполняем векторы данными для отправки
// right_upper_send[0] = Ex(Nx - 1, Ny - 1); // Ex
// right_upper_send[1] = Ey(Nx - 1, Ny - 1); // Ey
// right_upper_send[2] = Ez(Nx - 1, Ny - 1); // Ez
// right_upper_send[3] = Bx(Nx - 1, Ny - 1); // Bx
// right_upper_send[4] = By(Nx - 1, Ny - 1); // By
// right_upper_send[5] = Bz(Nx - 1, Ny - 1); // Bz

// int right_upper_coords[2] = {(coords[0] + 1) % dims[0], (coords[1] + 1) % dims[1]};
// int left_lower_coords[2] = {(coords[0] - 1 + dims[0]) % dims[0], (coords[1] - 1 + dims[1]) % dims[1]};
// int right_upper_rank;
// int left_lower_rank;
// MPI_Cart_rank(cart_comm, right_upper_coords, &right_upper_rank);
// MPI_Cart_rank(cart_comm, left_lower_coords, &left_lower_rank);

// // Отправляем и получаем данные
// MPI_Sendrecv(right_upper_send.data(), 6, MPI_DOUBLE, right_upper_rank, 6,
//  left_lower_receive.data(), 6, MPI_DOUBLE, left_lower_rank, 6, cart_comm, MPI_STATUS_IGNORE);

// // Записываем данные из вектора приёма в поля
// Ex(-1, -1) = left_lower_receive[0]; // Ex
// Ey(-1, -1) = left_lower_receive[1]; // Ey
// Ez(-1, -1) = left_lower_receive[2]; // Ez
// Bx(-1, -1) = left_lower_receive[3]; // Bx
// By(-1, -1) = left_lower_receive[4]; // By
// Bz(-1, -1) = left_lower_receive[5]; // Bz

// // ПУНКТ 7 ---------------------------------------------------

// // Векторы для отправки и приёма угловых элементов
// std::vector<double> right_lower_send(6, 0.0);
// std::vector<double> left_upper_receive(6, 0.0);

// // Заполняем векторы данными для отправки
// right_lower_send[0] = Ex(Nx - 1, 0); // Ex
// right_lower_send[1] = Ey(Nx - 1, 0); // Ey
// right_lower_send[2] = Ez(Nx - 1, 0); // Ez
// right_lower_send[3] = Bx(Nx - 1, 0); // Bx
// right_lower_send[4] = By(Nx - 1, 0); // By
// right_lower_send[5] = Bz(Nx - 1, 0); // Bz

// // Отправляем и получаем данные
// MPI_Sendrecv(right_lower_send.data(), 6, MPI_DOUBLE, right_lower_rank, 7,
//  left_upper_receive.data(), 6, MPI_DOUBLE, left_upper_rank, 7, cart_comm, MPI_STATUS_IGNORE);

// // Записываем данные из вектора приёма в поля
// Ex(-1, Ny) = left_upper_receive[0]; // Ex
// Ey(-1, Ny) = left_upper_receive[1]; // Ey
// Ez(-1, Ny) = left_upper_receive[2]; // Ez
// Bx(-1, Ny) = left_upper_receive[3]; // Bx
// By(-1, Ny) = left_upper_receive[4]; // By
// Bz(-1, Ny) = left_upper_receive[5]; // Bz

// // ПУНКТ 8 ---------------------------------------------------

// // Векторы для отправки и приёма угловых элементов
// std::vector<double> left_lower_send(6, 0.0);
// std::vector<double> right_upper_receive(6, 0.0);

// // Заполняем векторы данными для отправки
// left_lower_send[0] = Ex(0, 0); // Ex
// left_lower_send[1] = Ey(0, 0); // Ey
// left_lower_send[2] = Ez(0, 0); // Ez
// left_lower_send[3] = Bx(0, 0); // Bx
// left_lower_send[4] = By(0, 0); // By
// left_lower_send[5] = Bz(0, 0); // Bz

// // Отправляем и получаем данные
// MPI_Sendrecv(left_lower_send.data(), 6, MPI_DOUBLE, left_lower_rank, 8,
//  right_upper_receive.data(), 6, MPI_DOUBLE, right_upper_rank, 8, cart_comm, MPI_STATUS_IGNORE);

// // Записываем данные из вектора приёма в поля
// Ex(Nx, Ny) = right_upper_receive[0]; // Ex
// Ey(Nx, Ny) = right_upper_receive[1]; // Ey
// Ez(Nx, Ny) = right_upper_receive[2]; // Ez
// Bx(Nx, Ny) = right_upper_receive[3]; // Bx
// By(Nx, Ny) = right_upper_receive[4]; // By
// Bz(Nx, Ny) = right_upper_receive[5]; // Bz





// ---------------------------------------------------------------------------------------------------------------------
// Это для разбиения по одной оси OY (по строкам)

  // // Сначала рассматривается разбиение по OY для 2D

  // // Слева и справа (граничные условия для 2D) (Всё это внутри процесса, поэтому без MPI сообщений)
  // for (int64_t i = 0; i < Ny; ++i) {
  //   // Слева
  //   Ex(-1, i) = Ex(Nx - 1, i);
  //   Ey(-1, i) = Ey(Nx - 1, i);
  //   Ez(-1, i) = Ez(Nx - 1, i);
  //   Bx(-1, i) = Bx(Nx - 1, i);
  //   By(-1, i) = By(Nx - 1, i);
  //   Bz(-1, i) = Bz(Nx - 1, i);

  //   // Справа
  //   Ex(Nx, i) = Ex(0, i);
  //   Ey(Nx, i) = Ey(0, i);
  //   Ez(Nx, i) = Ez(0, i);
  //   Bx(Nx, i) = Bx(0, i);
  //   By(Nx, i) = By(0, i);
  //   Bz(Nx, i) = Bz(0, i);
  // }

  // // В NY фактически хранится Ny_local (то есть размер без оболочки)

  // // Нужно потом будет учесть ситуацию, когда только один процесс
  // // ТАКЖЕ ИНДЕКСЫ ЗДЕСЬ НЕ УЧИТЫВАЮТ OZ НАПРАВЛЕНИЕ!
  // // НЕ забыть добавить правильную обработку угловых значений для 0 и последнего процессов
  // // ****************************************************************************************************
  // // Низ каждого процесса передаётся вверх предыдущему процессу (P.S. не нулевая строка, а первая)

  // if (mpi_comm_size > 1) {

  //   if (rank != 0) {
  //     MPI_Send(&Ex(-1, 0), Nx + 2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
  //     MPI_Send(&Ey(-1, 0), Nx + 2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
  //     MPI_Send(&Ez(-1, 0), Nx + 2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
  //     MPI_Send(&Bx(-1, 0), Nx + 2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
  //     MPI_Send(&By(-1, 0), Nx + 2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
  //     MPI_Send(&Bz(-1, 0), Nx + 2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
  //   }

  //   // Принятие верха от низа следующего процесса (P.S. последняя строка)
  //   if (rank != mpi_comm_size - 1) {
  //     MPI_Recv(&Ex(-1, Ny), Nx + 2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Recv(&Ey(-1, Ny), Nx + 2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Recv(&Ez(-1, Ny), Nx + 2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Recv(&Bx(-1, Ny), Nx + 2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Recv(&By(-1, Ny), Nx + 2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Recv(&Bz(-1, Ny), Nx + 2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //   }
  //   // ****************************************************************************************************

  //   // Верх каждого процесса передаётся вниз следующему процессу (P.S. не последняя строка, а предпоследняя)
  //   if (rank != mpi_comm_size - 1) {
  //     MPI_Send(&Ex(-1, Ny - 1), Nx + 2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
  //     MPI_Send(&Ey(-1, Ny - 1), Nx + 2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
  //     MPI_Send(&Ez(-1, Ny - 1), Nx + 2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
  //     MPI_Send(&Bx(-1, Ny - 1), Nx + 2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
  //     MPI_Send(&By(-1, Ny - 1), Nx + 2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
  //     MPI_Send(&Bz(-1, Ny - 1), Nx + 2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
  //   }
    
  //   // Принятие низа от верха предыдущего процесса (P.S нулевая строка)
  //   if (rank != 0) {
  //     MPI_Recv(Ex.data(), Nx + 2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Recv(Ey.data(), Nx + 2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Recv(Ez.data(), Nx + 2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Recv(Bx.data(), Nx + 2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Recv(By.data(), Nx + 2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Recv(Bz.data(), Nx + 2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //   }

  //   // ****************************************************************************************************

  //   // Синхронизация для нулевого и последнего процессов

  //   // Отправка низа вверх последнего процесса (P.S. не нулевая строка, а первая), а также получение низа (нулевой строки) от последнего процесса
  //   if (rank == 0) {
  //     MPI_Sendrecv(&Ex(-1, 0), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 0, Ex.data(), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Sendrecv(&Ey(-1, 0), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 0, Ey.data(), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Sendrecv(&Ez(-1, 0), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 0, Ez.data(), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Sendrecv(&Bx(-1, 0), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 0, Bx.data(), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Sendrecv(&By(-1, 0), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 0, By.data(), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Sendrecv(&Bz(-1, 0), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 0, Bz.data(), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //   }

  //   // Отправка верха вниз нулевого процесса (P.S. не последняя строка, а предпоследняя), а также получение верха (последней строки) от нулевого процесса
  //   if (rank == mpi_comm_size - 1) {
  //     MPI_Sendrecv(&Ex(-1, Ny - 1), Nx + 2, MPI_DOUBLE, 0, 1, &Ex(-1, Ny), Nx + 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Sendrecv(&Ey(-1, Ny - 1), Nx + 2, MPI_DOUBLE, 0, 1, &Ey(-1, Ny), Nx + 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Sendrecv(&Ez(-1, Ny - 1), Nx + 2, MPI_DOUBLE, 0, 1, &Ez(-1, Ny), Nx + 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Sendrecv(&Bx(-1, Ny - 1), Nx + 2, MPI_DOUBLE, 0, 1, &Bx(-1, Ny), Nx + 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Sendrecv(&By(-1, Ny - 1), Nx + 2, MPI_DOUBLE, 0, 1, &By(-1, Ny), Nx + 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Sendrecv(&Bz(-1, Ny - 1), Nx + 2, MPI_DOUBLE, 0, 1, &Bz(-1, Ny), Nx + 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //   }

  // } else {
  //   // Верх (с угловыми узлами)
  //   for (int64_t x = -1; x < Nx + 1; ++x) {
  //     Ex(x, Ny) = Ex(x, 0);
  //     Ey(x, Ny) = Ey(x, 0);
  //     Ez(x, Ny) = Ez(x, 0);
  //     Bx(x, Ny) = Bx(x, 0);
  //     By(x, Ny) = By(x, 0);
  //     Bz(x, Ny) = Bz(x, 0);
  //   }
  //   // Низ (с угловыми узлами)
  //   for (int64_t x = -1; x < Nx + 1; ++x) {
  //     Ex(x, -1) = Ex(x, Ny - 1);
  //     Ey(x, -1) = Ey(x, Ny - 1);
  //     Ez(x, -1) = Ez(x, Ny - 1);
  //     Bx(x, -1) = Bx(x, Ny - 1);
  //     By(x, -1) = By(x, Ny - 1);
  //     Bz(x, -1) = Bz(x, Ny - 1);
  //   }
  // }

// ---------------------------------------------------------------------------------------------------------------------


  // ****************************************************************************************************



//   // Синхронизация границ для 2д случая
//   // Угловые узлы обрабатываются сразу, а не отдельно

//   // Слева (за исключением угловых узлов)
// #pragma omp for
//   for (int64_t y = 0; y < Ny; ++y) {
//     Ex(-1, y) = Ex(Nx - 1, y);
//     Ey(-1, y) = Ey(Nx - 1, y);
//     Ez(-1, y) = Ez(Nx - 1, y);
//     Bx(-1, y) = Bx(Nx - 1, y);
//     By(-1, y) = By(Nx - 1, y);
//     Bz(-1, y) = Bz(Nx - 1, y);
//   }
//   // Справа (за исключением угловых узлов)
// #pragma omp for
//   for (int64_t y = 0; y < Ny; ++y) {
//     Ex(Nx, y) = Ex(0, y);
//     Ey(Nx, y) = Ey(0, y);
//     Ez(Nx, y) = Ez(0, y);
//     Bx(Nx, y) = Bx(0, y);
//     By(Nx, y) = By(0, y);
//     Bz(Nx, y) = Bz(0, y);
//   }
// #pragma omp for
//   for (int64_t x = -1; x < Nx + 1; ++x) {
//     Ex(x, Ny) = Ex(x, 0);
//     Ey(x, Ny) = Ey(x, 0);
//     Ez(x, Ny) = Ez(x, 0);
//     Bx(x, Ny) = Bx(x, 0);
//     By(x, Ny) = By(x, 0);
//     Bz(x, Ny) = Bz(x, 0);
//   }
//   // Низ (с угловыми узлами)
// #pragma omp for
//   for (int64_t x = -1; x < Nx + 1; ++x) {
//     Ex(x, -1) = Ex(x, Ny - 1);
//     Ey(x, -1) = Ey(x, Ny - 1);
//     Ez(x, -1) = Ez(x, Ny - 1);
//     Bx(x, -1) = Bx(x, Ny - 1);
//     By(x, -1) = By(x, Ny - 1);
//     Bz(x, -1) = Bz(x, Ny - 1);
//   }
}

void FDTD::FDTD::boundary_synchronization_3D(MPI_Comm cart_comm) {

  int rank = 0;
  int mpi_comm_size = 0;
  MPI_Comm_rank(cart_comm, &rank);
  MPI_Comm_size(cart_comm, &mpi_comm_size);

  // Важно:
  // Может быть ситуация, когда процессов для 3D не будет хватать
  // Пример:
    // MPI_dimension (3D or 2D domain): 3D
    // Count of processes: 10
    // Ox processes: 2 Oy processes: 5 Oz processes: 1
  // Но вроде ошибки не должно быть в моём коде, если такой случай. Пока хз, нужно будет в конечном счёте это всё тестировать.

  // Инструкция по синхронизации границ для 3D случая в случае разбиения по всем трём осям:
  // Угловые граничный элементы будут учитываться сразу, а не отдельно
  // ---------------------------------------------------------------------------------------------------------------------
  // ЧИТАТЬ СЛЕДУЮЩУЮ ИНФОРМАЦИЮ, ПРЕДСТАВЛЯЯ ВИД СБОКУ НА КУБИКИ
  // ---------------------------------------------------------------------------------------------------------------------
  
  // ТЕГ 1. Сначала отправляется левая грань каждого процесса вправо левому процессу.
  //    и принятие левой грани от правого процесса
  // ТЕГ 2. Отправка правой грани каждого процесса влево правому процессу
  //    и принятие правой грани от левого процесса
  // ТЕГ 3. Отправка дальней грани каждого процесса на ближнюю грань дальнего (следующего) процесса
  //    и принятие дальней грани от ближнего процесса
  // ТЕГ 4. Отправка ближней грани каждого процесса на дальнюю грань ближнего (следующего) процесса
  //    и принятие ближней грани от дальнего процесса
  // ТЕГ 5. Отправка верхней грани каждого процесса вниз верхнему процессу
  //    и принятие верхней грани от нижнего процесса
  // ТЕГ 6. Отправка нижней грани каждого процесса вверх нижнему процессу
  //    и принятие нижней грани от верхнего процесса

  // ВИДИМО, для того, чтобы включить обработку угловых узлов, для это нужно включить узлы противоположных сторон (левая-правая или дальняя-ближняя)

  int left, right, up, down, front, back;
  // MPI_Cart_shift(cart_comm, 0, 1, &up, &down);
  // MPI_Cart_shift(cart_comm, 1, 1, &left, &right);
  // MPI_Cart_shift(cart_comm, 2, 1, &back, &front);

MPI_Cart_shift(cart_comm, 0, 1, &left, &right); 
MPI_Cart_shift(cart_comm, 1, 1, &back, &front);
MPI_Cart_shift(cart_comm, 2, 1, &down, &up);

  // std::cout << std::endl;
  // std::cout << "Rank: " << rank << " Left: " << left << " Right: " << right << " Up: " << up << " Down: " << down << " Front: " << front << " Back: " << back 
  //  << std::endl;

  // 1 ПУНКТ ---------------------------------------------------
  // Векторы для отправки и приёма левой грани
  std::vector<std::vector<double>> left_send(6, std::vector<double>(Ny * Nz, 0.0));
  std::vector<std::vector<double>> right_receive(6, std::vector<double>(Ny * Nz, 0.0));

  // Заполняем векторы данными для отправки
  for (int64_t y = 0; y < Ny; ++y) {
    for (int64_t z = 0; z < Nz; ++z) {
      int index = y * Nz + z;
      left_send[0][index] = Ex(0, y, z); // Ex
      left_send[1][index] = Ey(0, y, z); // Ey
      left_send[2][index] = Ez(0, y, z); // Ez
      left_send[3][index] = Bx(0, y, z); // Bx
      left_send[4][index] = By(0, y, z); // By
      left_send[5][index] = Bz(0, y, z); // Bz
    }
  }

  // Отправляем и получаем данные
  MPI_Sendrecv(left_send[0].data(), Ny * Nz, MPI_DOUBLE, left, 1,
    right_receive[0].data(), Ny * Nz, MPI_DOUBLE, right, 1, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(left_send[1].data(), Ny * Nz, MPI_DOUBLE, left, 1,
    right_receive[1].data(), Ny * Nz, MPI_DOUBLE, right, 1, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(left_send[2].data(), Ny * Nz, MPI_DOUBLE, left, 1,
    right_receive[2].data(), Ny * Nz, MPI_DOUBLE, right, 1, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(left_send[3].data(), Ny * Nz, MPI_DOUBLE, left, 1,
    right_receive[3].data(), Ny * Nz, MPI_DOUBLE, right, 1, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(left_send[4].data(), Ny * Nz, MPI_DOUBLE, left, 1,
    right_receive[4].data(), Ny * Nz, MPI_DOUBLE, right, 1, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(left_send[5].data(), Ny * Nz, MPI_DOUBLE, left, 1,
    right_receive[5].data(), Ny * Nz, MPI_DOUBLE, right, 1, cart_comm, MPI_STATUS_IGNORE);

  // Записываем данные из вектора приёма в поля
  for (int64_t y = 0; y < Ny; ++y) {
    for (int64_t z = 0; z < Nz; ++z) {
      Ex(Nx, y, z) = right_receive[0][y * Nz + z]; // Ex
      Ey(Nx, y, z) = right_receive[1][y * Nz + z]; // Ey
      Ez(Nx, y, z) = right_receive[2][y * Nz + z]; // Ez
      Bx(Nx, y, z) = right_receive[3][y * Nz + z]; // Bx
      By(Nx, y, z) = right_receive[4][y * Nz + z]; // By
      Bz(Nx, y, z) = right_receive[5][y * Nz + z]; // Bz
    }
  }
  
  // 2 ПУНКТ ---------------------------------------------------
  // Векторы для отправки и приёма правой грани
  std::vector<std::vector<double>> right_send(6, std::vector<double>(Ny * Nz, 0.0));
  std::vector<std::vector<double>> left_receive(6, std::vector<double>(Ny * Nz, 0.0));

  // Заполняем векторы данными для отправки
  for (int64_t y = 0; y < Ny; ++y) {
    for (int64_t z = 0; z < Nz; ++z) {
      int index = y * Nz + z;
      right_send[0][index] = Ex(Nx - 1, y, z); // Ex
      right_send[1][index] = Ey(Nx - 1, y, z); // Ey
      right_send[2][index] = Ez(Nx - 1, y, z); // Ez
      right_send[3][index] = Bx(Nx - 1, y, z); // Bx
      right_send[4][index] = By(Nx - 1, y, z); // By
      right_send[5][index] = Bz(Nx - 1, y, z); // Bz
    }
  }

  // Отправляем и получаем данные
  MPI_Sendrecv(right_send[0].data(), Ny * Nz, MPI_DOUBLE, right, 2,
    left_receive[0].data(), Ny * Nz, MPI_DOUBLE, left, 2, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(right_send[1].data(), Ny * Nz, MPI_DOUBLE, right, 2,
    left_receive[1].data(), Ny * Nz, MPI_DOUBLE, left, 2, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(right_send[2].data(), Ny * Nz, MPI_DOUBLE, right, 2,
    left_receive[2].data(), Ny * Nz, MPI_DOUBLE, left, 2, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(right_send[3].data(), Ny * Nz, MPI_DOUBLE, right, 2,
    left_receive[3].data(), Ny * Nz, MPI_DOUBLE, left, 2, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(right_send[4].data(), Ny * Nz, MPI_DOUBLE, right, 2,
    left_receive[4].data(), Ny * Nz, MPI_DOUBLE, left, 2, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(right_send[5].data(), Ny * Nz, MPI_DOUBLE, right, 2,
    left_receive[5].data(), Ny * Nz, MPI_DOUBLE, left, 2, cart_comm, MPI_STATUS_IGNORE);

  // Записываем данные из вектора приёма в поля
  for (int64_t y = 0; y < Ny; ++y) {
    for (int64_t z = 0; z < Nz; ++z) {
      Ex(-1, y, z) = left_receive[0][y * Nz + z]; // Ex
      Ey(-1, y, z) = left_receive[1][y * Nz + z]; // Ey
      Ez(-1, y, z) = left_receive[2][y * Nz + z]; // Ez
      Bx(-1, y, z) = left_receive[3][y * Nz + z]; // Bx
      By(-1, y, z) = left_receive[4][y * Nz + z]; // By
      Bz(-1, y, z) = left_receive[5][y * Nz + z]; // Bz
    }
  }
  
  // 3 ПУНКТ ---------------------------------------------------
  // Векторы для отправки и приёма дальней грани
  std::vector<std::vector<double>> back_send(6, std::vector<double>((Nx + 2) * Nz, 0.0));
  std::vector<std::vector<double>> front_receive(6, std::vector<double>((Nx + 2) * Nz, 0.0));

  // Заполняем векторы данными для отправки
  for (int64_t x = -1; x < Nx + 1; ++x) {
    for (int64_t z = 0; z < Nz; ++z) {
      int index = (x + 1) * Nz + z;
      back_send[0][index] = Ex(x, Ny - 1, z); // Ex
      back_send[1][index] = Ey(x, Ny - 1, z); // Ey
      back_send[2][index] = Ez(x, Ny - 1, z); // Ez
      back_send[3][index] = Bx(x, Ny - 1, z); // Bx
      back_send[4][index] = By(x, Ny - 1, z); // By
      back_send[5][index] = Bz(x, Ny - 1, z); // Bz
    }
  }

  // Отправляем и получаем данные
  MPI_Sendrecv(back_send[0].data(), (Nx + 2) * Nz, MPI_DOUBLE, back, 3,
    front_receive[0].data(), (Nx + 2) * Nz, MPI_DOUBLE, front, 3, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(back_send[1].data(), (Nx + 2) * Nz, MPI_DOUBLE, back, 3,
    front_receive[1].data(), (Nx + 2) * Nz, MPI_DOUBLE, front, 3, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(back_send[2].data(), (Nx + 2) * Nz, MPI_DOUBLE, back, 3,
    front_receive[2].data(), (Nx + 2) * Nz, MPI_DOUBLE, front, 3, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(back_send[3].data(), (Nx + 2) * Nz, MPI_DOUBLE, back, 3,
    front_receive[3].data(), (Nx + 2) * Nz, MPI_DOUBLE, front, 3, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(back_send[4].data(), (Nx + 2) * Nz, MPI_DOUBLE, back, 3,
    front_receive[4].data(), (Nx + 2) * Nz, MPI_DOUBLE, front, 3, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(back_send[5].data(), (Nx + 2) * Nz, MPI_DOUBLE, back, 3,
    front_receive[5].data(), (Nx + 2) * Nz, MPI_DOUBLE, front, 3, cart_comm, MPI_STATUS_IGNORE);

  // Записываем данные из вектора приёма в поля
  for (int64_t x = -1; x < Nx + 1; ++x) {
    for (int64_t z = 0; z < Nz; ++z) {
      Ex(x, -1, z) = front_receive[0][(x + 1) * Nz + z]; // Ex
      Ey(x, -1, z) = front_receive[1][(x + 1) * Nz + z]; // Ey
      Ez(x, -1, z) = front_receive[2][(x + 1) * Nz + z]; // Ez
      Bx(x, -1, z) = front_receive[3][(x + 1) * Nz + z]; // Bx
      By(x, -1, z) = front_receive[4][(x + 1) * Nz + z]; // By
      Bz(x, -1, z) = front_receive[5][(x + 1) * Nz + z]; // Bz
    }
  }

  // 4 ПУНКТ ---------------------------------------------------
  // Векторы для отправки и приёма ближней грани
  std::vector<std::vector<double>> front_send(6, std::vector<double>((Nx + 2) * Nz, 0.0));
  std::vector<std::vector<double>> back_receive(6, std::vector<double>((Nx + 2) * Nz, 0.0));

  // Заполняем векторы данными для отправки
  for (int64_t x = -1; x < Nx + 1; ++x) {
    for (int64_t z = 0; z < Nz; ++z) {
      int index = (x + 1) * Nz + z;
      front_send[0][index] = Ex(x, 0, z); // Ex
      front_send[1][index] = Ey(x, 0, z); // Ey
      front_send[2][index] = Ez(x, 0, z); // Ez
      front_send[3][index] = Bx(x, 0, z); // Bx
      front_send[4][index] = By(x, 0, z); // By
      front_send[5][index] = Bz(x, 0, z); // Bz
    }
  }

  // Отправляем и получаем данные
  MPI_Sendrecv(front_send[0].data(), (Nx + 2) * Nz, MPI_DOUBLE, front, 4,
    back_receive[0].data(), (Nx + 2) * Nz, MPI_DOUBLE, back, 4, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(front_send[1].data(), (Nx + 2) * Nz, MPI_DOUBLE, front, 4,
    back_receive[1].data(), (Nx + 2) * Nz, MPI_DOUBLE, back, 4, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(front_send[2].data(), (Nx + 2) * Nz, MPI_DOUBLE, front, 4,
    back_receive[2].data(), (Nx + 2) * Nz, MPI_DOUBLE, back, 4, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(front_send[3].data(), (Nx + 2) * Nz, MPI_DOUBLE, front, 4,
    back_receive[3].data(), (Nx + 2) * Nz, MPI_DOUBLE, back, 4, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(front_send[4].data(), (Nx + 2) * Nz, MPI_DOUBLE, front, 4,
    back_receive[4].data(), (Nx + 2) * Nz, MPI_DOUBLE, back, 4, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(front_send[5].data(), (Nx + 2) * Nz, MPI_DOUBLE, front, 4,
    back_receive[5].data(), (Nx + 2) * Nz, MPI_DOUBLE, back, 4, cart_comm, MPI_STATUS_IGNORE);

  // Записываем данные из вектора приёма в поля
  for (int64_t x = -1; x < Nx + 1; ++x) {
    for (int64_t z = 0; z < Nz; ++z) {
      Ex(x, Ny, z) = back_receive[0][(x + 1) * Nz + z]; // Ex
      Ey(x, Ny, z) = back_receive[1][(x + 1) * Nz + z]; // Ey
      Ez(x, Ny, z) = back_receive[2][(x + 1) * Nz + z]; // Ez
      Bx(x, Ny, z) = back_receive[3][(x + 1) * Nz + z]; // Bx
      By(x, Ny, z) = back_receive[4][(x + 1) * Nz + z]; // By
      Bz(x, Ny, z) = back_receive[5][(x + 1) * Nz + z]; // Bz
    }
  }

  // 5 ПУНКТ ---------------------------------------------------
  // Векторы для отправки и приёма верхней грани
  std::vector<std::vector<double>> up_send(6, std::vector<double>((Nx + 2) * (Ny + 2), 0.0));
  std::vector<std::vector<double>> down_receive(6, std::vector<double>((Nx + 2) * (Ny + 2), 0.0));

  // Заполняем векторы данными для отправки
  for (int64_t x = -1; x < Nx + 1; ++x) {
    for (int64_t y = -1; y < Ny + 1; ++y) {
      int index = (x + 1) * (Ny + 2) + y + 1;
      up_send[0][index] = Ex(x, y, Nz - 1); // Ex
      up_send[1][index] = Ey(x, y, Nz - 1); // Ey
      up_send[2][index] = Ez(x, y, Nz - 1); // Ez
      up_send[3][index] = Bx(x, y, Nz - 1); // Bx
      up_send[4][index] = By(x, y, Nz - 1); // By
      up_send[5][index] = Bz(x, y, Nz - 1); // Bz
    }
  }

  // Отправляем и получаем данные
  MPI_Sendrecv(up_send[0].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 5,
    down_receive[0].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 5, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(up_send[1].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 5,
    down_receive[1].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 5, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(up_send[2].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 5,
    down_receive[2].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 5, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(up_send[3].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 5,
    down_receive[3].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 5, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(up_send[4].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 5,
    down_receive[4].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 5, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(up_send[5].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 5,
    down_receive[5].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 5, cart_comm, MPI_STATUS_IGNORE);

  // Записываем данные из вектора приёма в поля
  for (int64_t x = -1; x < Nx + 1; ++x) {
    for (int64_t y = -1; y < Ny + 1; ++y) {
      Ex(x, y, -1) = down_receive[0][(x + 1) * (Ny + 2) + y + 1]; // Ex
      Ey(x, y, -1) = down_receive[1][(x + 1) * (Ny + 2) + y + 1]; // Ey
      Ez(x, y, -1) = down_receive[2][(x + 1) * (Ny + 2) + y + 1]; // Ez
      Bx(x, y, -1) = down_receive[3][(x + 1) * (Ny + 2) + y + 1]; // Bx
      By(x, y, -1) = down_receive[4][(x + 1) * (Ny + 2) + y + 1]; // By
      Bz(x, y, -1) = down_receive[5][(x + 1) * (Ny + 2) + y + 1]; // Bz
    }
  }

  // 6 ПУНКТ ---------------------------------------------------
  // Векторы для отправки и приёма нижней грани
  std::vector<std::vector<double>> down_send(6, std::vector<double>((Nx + 2) * (Ny + 2), 0.0));
  std::vector<std::vector<double>> up_receive(6, std::vector<double>((Nx + 2) * (Ny + 2), 0.0));

  // Заполняем векторы данными для отправки
  for (int64_t x = -1; x < Nx + 1; ++x) {
    for (int64_t y = -1; y < Ny + 1; ++y) {
      int index = (x + 1) * (Ny + 2) + y + 1;
      down_send[0][index] = Ex(x, y, 0); // Ex
      down_send[1][index] = Ey(x, y, 0); // Ey
      down_send[2][index] = Ez(x, y, 0); // Ez
      down_send[3][index] = Bx(x, y, 0); // Bx
      down_send[4][index] = By(x, y, 0); // By
      down_send[5][index] = Bz(x, y, 0); // Bz
    }
  }

  // Отправляем и получаем данные
  MPI_Sendrecv(down_send[0].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 6,
    up_receive[0].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 6, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(down_send[1].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 6,
    up_receive[1].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 6, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(down_send[2].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 6,
    up_receive[2].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 6, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(down_send[3].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 6,
    up_receive[3].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 6, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(down_send[4].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 6,
    up_receive[4].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 6, cart_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(down_send[5].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 6,
    up_receive[5].data(), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 6, cart_comm, MPI_STATUS_IGNORE);

  // Записываем данные из вектора приёма в поля
  for (int64_t x = -1; x < Nx + 1; ++x) {
    for (int64_t y = -1; y < Ny + 1; ++y) {
      Ex(x, y, Nz) = up_receive[0][(x + 1) * (Ny + 2) + y + 1]; // Ex
      Ey(x, y, Nz) = up_receive[1][(x + 1) * (Ny + 2) + y + 1]; // Ey
      Ez(x, y, Nz) = up_receive[2][(x + 1) * (Ny + 2) + y + 1]; // Ez
      Bx(x, y, Nz) = up_receive[3][(x + 1) * (Ny + 2) + y + 1]; // Bx
      By(x, y, Nz) = up_receive[4][(x + 1) * (Ny + 2) + y + 1]; // By
      Bz(x, y, Nz) = up_receive[5][(x + 1) * (Ny + 2) + y + 1]; // Bz
    }
  }

  


// ---------------------------------------------------------------------------------------------------------------------
// ---------------------Не рабочий вариант ---------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

  // // 1 ПУНКТ ---------------------------------------------------
  // // Векторы для отправки и приёма левой боковой плоскости
  // std::vector<std::vector<double>> left_send(6, std::vector<double>((Ny + 2) * (Nz + 2), 0.0));
  // std::vector<std::vector<double>> right_receive(6, std::vector<double>((Ny + 2) * (Nz + 2), 0.0));

  // // Заполняем векторы данными для отправки
  // for (int64_t y = 0; y < Ny + 2; ++y) {
  //   for (int64_t z = 0; z < Nz + 2; ++z) {
  //     int index = y * (Nz + 2) + z;
  //     left_send[0][index] = Ex(0, y - 1, z - 1); // Ex
  //     left_send[1][index] = Ey(0, y - 1, z - 1); // Ey
  //     left_send[2][index] = Ez(0, y - 1, z - 1); // Ez
  //     left_send[3][index] = Bx(0, y - 1, z - 1); // Bx
  //     left_send[4][index] = By(0, y - 1, z - 1); // By
  //     left_send[5][index] = Bz(0, y - 1, z - 1); // Bz
  //   }
  // }

  // // Отправляем и получаем данные
  // MPI_Sendrecv(left_send[0].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, left, 1,
  //   right_receive[0].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, right, 1, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(left_send[1].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, left, 1, 
  //   right_receive[1].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, right, 1, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(left_send[2].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, left, 1, 
  //   right_receive[2].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, right, 1, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(left_send[3].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, left, 1,
  //   right_receive[3].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, right, 1, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(left_send[4].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, left, 1,
  //   right_receive[4].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, right, 1, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(left_send[5].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, left, 1,
  //   right_receive[5].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, right, 1, cart_comm, MPI_STATUS_IGNORE);

  // // Записываем данные из вектора приёма в поля
  // for (int64_t y = 0; y < Ny + 2; ++y) {
  //   for (int64_t z = 0; z < Nz + 2; ++z) {
  //     Ex(Nx, y - 1, z - 1) = right_receive[0][y * Nz + z]; // Ex
  //     Ey(Nx, y - 1, z - 1) = right_receive[1][y * Nz + z]; // Ey
  //     Ez(Nx, y - 1, z - 1) = right_receive[2][y * Nz + z]; // Ez
  //     Bx(Nx, y - 1, z - 1) = right_receive[3][y * Nz + z]; // Bx
  //     By(Nx, y - 1, z - 1) = right_receive[4][y * Nz + z]; // By
  //     Bz(Nx, y - 1, z - 1) = right_receive[5][y * Nz + z]; // Bz
  //   }
  // }

  // // 2 ПУНКТ ---------------------------------------------------
  // // Векторы для отправки и приёма правой боковой плоскости
  // std::vector<std::vector<double>> right_send(6, std::vector<double>((Ny + 2) * (Nz + 2), 0.0));
  // std::vector<std::vector<double>> left_receive(6, std::vector<double>((Ny + 2) * (Nz + 2), 0.0));

  // // Заполняем векторы данными для отправки
  // for (int64_t y = 0; y < Ny + 2; ++y) {
  //   for (int64_t z = 0; z < Nz + 2; ++z) {
  //     right_send[0][y * Nz + z] = Ex(Nx - 1, y - 1, z - 1); // Ex
  //     right_send[1][y * Nz + z] = Ey(Nx - 1, y - 1, z - 1); // Ey
  //     right_send[2][y * Nz + z] = Ez(Nx - 1, y - 1, z - 1); // Ez
  //     right_send[3][y * Nz + z] = Bx(Nx - 1, y - 1, z - 1); // Bx
  //     right_send[4][y * Nz + z] = By(Nx - 1, y - 1, z - 1); // By
  //     right_send[5][y * Nz + z] = Bz(Nx - 1, y - 1, z - 1); // Bz
  //   }
  // }

  // // Отправляем и получаем данные
  // MPI_Sendrecv(right_send[0].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, right, 2,
  //   left_receive[0].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, left, 2, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(right_send[1].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, right, 2,
  //   left_receive[1].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, left, 2, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(right_send[2].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, right, 2,
  //   left_receive[2].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, left, 2, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(right_send[3].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, right, 2,
  //   left_receive[3].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, left, 2, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(right_send[4].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, right, 2,
  //   left_receive[4].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, left, 2, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(right_send[5].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, right, 2,
  //   left_receive[5].data(), (Ny + 2) * (Nz + 2), MPI_DOUBLE, left, 2, cart_comm, MPI_STATUS_IGNORE);

  // // Записываем данные из вектора приёма в поля
  // for (int64_t y = 0; y < Ny + 2; ++y) {
  //   for (int64_t z = 0; z < Nz + 2; ++z) {
  //     Ex(-1, y - 1, z - 1) = left_receive[0][y * Nz + z]; // Ex
  //     Ey(-1, y - 1, z - 1) = left_receive[1][y * Nz + z]; // Ey
  //     Ez(-1, y - 1, z - 1) = left_receive[2][y * Nz + z]; // Ez
  //     Bx(-1, y - 1, z - 1) = left_receive[3][y * Nz + z]; // Bx
  //     By(-1, y - 1, z - 1) = left_receive[4][y * Nz + z]; // By
  //     Bz(-1, y - 1, z - 1) = left_receive[5][y * Nz + z]; // Bz
  //   }
  // }

  // // 3 ПУНКТ ---------------------------------------------------
  // // Векторы для отправки и приёма дальней плоскости
  // std::vector<std::vector<double>> back_send(6, std::vector<double>(Nx * (Nz + 2), 0.0));
  // std::vector<std::vector<double>> front_receive(6, std::vector<double>(Nx * (Nz + 2), 0.0));

  // // Заполняем векторы данными для отправки
  // for (int64_t x = 0; x < Nx; ++x) {
  //   for (int64_t z = 0; z < Nz + 2; ++z) {
  //     back_send[0][x * Nz + z] = Ex(x, Ny - 1, z - 1); // Ex
  //     back_send[1][x * Nz + z] = Ey(x, Ny - 1, z - 1); // Ey
  //     back_send[2][x * Nz + z] = Ez(x, Ny - 1, z - 1); // Ez
  //     back_send[3][x * Nz + z] = Bx(x, Ny - 1, z - 1); // Bx
  //     back_send[4][x * Nz + z] = By(x, Ny - 1, z - 1); // By
  //     back_send[5][x * Nz + z] = Bz(x, Ny - 1, z - 1); // Bz
  //   }
  // }

  // // Отправляем и получаем данные
  // MPI_Sendrecv(back_send[0].data(), Nx * (Nz + 2), MPI_DOUBLE, back, 3,
  //   front_receive[0].data(), Nx * (Nz + 2), MPI_DOUBLE, front, 3, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(back_send[1].data(), Nx * (Nz + 2), MPI_DOUBLE, back, 3,
  //   front_receive[1].data(), Nx * (Nz + 2), MPI_DOUBLE, front, 3, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(back_send[2].data(), Nx * (Nz + 2), MPI_DOUBLE, back, 3,
  //   front_receive[2].data(), Nx * (Nz + 2), MPI_DOUBLE, front, 3, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(back_send[3].data(), Nx * (Nz + 2), MPI_DOUBLE, back, 3,
  //   front_receive[3].data(), Nx * (Nz + 2), MPI_DOUBLE, front, 3, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(back_send[4].data(), Nx * (Nz + 2), MPI_DOUBLE, back, 3,
  //   front_receive[4].data(), Nx * (Nz + 2), MPI_DOUBLE, front, 3, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(back_send[5].data(), Nx * (Nz + 2), MPI_DOUBLE, back, 3, 
  //   front_receive[5].data(), Nx * (Nz + 2), MPI_DOUBLE, front, 3, cart_comm, MPI_STATUS_IGNORE);

  // // Записываем данные из вектора приёма в поля
  // for (int64_t x = 0; x < Nx; ++x) {
  //   for (int64_t z = 0; z < Nz + 2; ++z) {
  //     Ex(x, -1, z - 1) = front_receive[0][x * Nz + z]; // Ex
  //     Ey(x, -1, z - 1) = front_receive[1][x * Nz + z]; // Ey
  //     Ez(x, -1, z - 1) = front_receive[2][x * Nz + z]; // Ez
  //     Bx(x, -1, z - 1) = front_receive[3][x * Nz + z]; // Bx
  //     By(x, -1, z - 1) = front_receive[4][x * Nz + z]; // By
  //     Bz(x, -1, z - 1) = front_receive[5][x * Nz + z]; // Bz
  //   }
  // }

  // // 4 ПУНКТ ---------------------------------------------------
  // // Векторы для отправки и приёма ближней плоскости
  // std::vector<std::vector<double>> front_send(6, std::vector<double>(Nx * (Nz + 2), 0.0));
  // std::vector<std::vector<double>> back_receive(6, std::vector<double>(Nx * (Nz + 2), 0.0));

  // // Заполняем векторы данными для отправки
  // for (int64_t x = 0; x < Nx; ++x) {
  //   for (int64_t z = 0; z < Nz + 2; ++z) {
  //     front_send[0][x * Nz + z] = Ex(x, 0, z - 1); // Ex
  //     front_send[1][x * Nz + z] = Ey(x, 0, z - 1); // Ey
  //     front_send[2][x * Nz + z] = Ez(x, 0, z - 1); // Ez
  //     front_send[3][x * Nz + z] = Bx(x, 0, z - 1); // Bx
  //     front_send[4][x * Nz + z] = By(x, 0, z - 1); // By
  //     front_send[5][x * Nz + z] = Bz(x, 0, z - 1); // Bz
  //   }
  // }

  // // Отправляем и получаем данные
  // MPI_Sendrecv(front_send[0].data(), Nx * (Nz + 2), MPI_DOUBLE, front, 4,
  //   back_receive[0].data(), Nx * (Nz + 2), MPI_DOUBLE, back, 4, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(front_send[1].data(), Nx * (Nz + 2), MPI_DOUBLE, front, 4,
  //   back_receive[1].data(), Nx * (Nz + 2), MPI_DOUBLE, back, 4, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(front_send[2].data(), Nx * (Nz + 2), MPI_DOUBLE, front, 4,
  //   back_receive[2].data(), Nx * (Nz + 2), MPI_DOUBLE, back, 4, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(front_send[3].data(), Nx * (Nz + 2), MPI_DOUBLE, front, 4,
  //   back_receive[3].data(), Nx * (Nz + 2), MPI_DOUBLE, back, 4, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(front_send[4].data(), Nx * (Nz + 2), MPI_DOUBLE, front, 4,
  //   back_receive[4].data(), Nx * (Nz + 2), MPI_DOUBLE, back, 4, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(front_send[5].data(), Nx * (Nz + 2), MPI_DOUBLE, front, 4,
  //   back_receive[5].data(), Nx * (Nz + 2), MPI_DOUBLE, back, 4, cart_comm, MPI_STATUS_IGNORE);
  
  // // Записываем данные из вектора приёма в поля
  // for (int64_t x = 0; x < Nx; ++x) {
  //   for (int64_t z = 0; z < Nz + 2; ++z) {
  //     Ex(x, Ny, z - 1) = back_receive[0][x * Nz + z]; // Ex
  //     Ey(x, Ny, z - 1) = back_receive[1][x * Nz + z]; // Ey
  //     Ez(x, Ny, z - 1) = back_receive[2][x * Nz + z]; // Ez
  //     Bx(x, Ny, z - 1) = back_receive[3][x * Nz + z]; // Bx
  //     By(x, Ny, z - 1) = back_receive[4][x * Nz + z]; // By
  //     Bz(x, Ny, z - 1) = back_receive[5][x * Nz + z]; // Bz
  //   }
  // }

  // // 5 ПУНКТ ---------------------------------------------------
  // // Отправляем и получаем данные
  // MPI_Sendrecv(&Ex(-1, -1, Nz - 1), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 5,
  //   &Ex(-1, -1, 0), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 5, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(&Ey(-1, -1, Nz - 1), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 5,
  //   &Ey(-1, -1, 0), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 5, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(&Ez(-1, -1, Nz - 1), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 5,
  //   &Ez(-1, -1, 0), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 5, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(&Bx(-1, -1, Nz - 1), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 5,
  //   &Bx(-1, -1, 0), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 5, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(&By(-1, -1, Nz - 1), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 5,
  //   &By(-1, -1, 0), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 5, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(&Bz(-1, -1, Nz - 1), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 5,
  //   &Bz(-1, -1, 0), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 5, cart_comm, MPI_STATUS_IGNORE);

  // // 6 ПУНКТ ---------------------------------------------------
  // // Отправляем и получаем данные
  // MPI_Sendrecv(&Ex(-1, -1, 0), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 6,
  //   &Ex(-1, -1, Nz), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 6, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(&Ey(-1, -1, 0), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 6,
  //   &Ey(-1, -1, Nz), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 6, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(&Ez(-1, -1, 0), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 6,
  //   &Ez(-1, -1, Nz), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 6, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(&Bx(-1, -1, 0), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 6,
  //   &Bx(-1, -1, Nz), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 6, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(&By(-1, -1, 0), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 6,
  //   &By(-1, -1, Nz), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 6, cart_comm, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(&Bz(-1, -1, 0), (Nx + 2) * (Ny + 2), MPI_DOUBLE, down, 6,
  //   &Bz(-1, -1, Nz), (Nx + 2) * (Ny + 2), MPI_DOUBLE, up, 6, cart_comm, MPI_STATUS_IGNORE);

// ---------------------------------------------------------------------------------------------------------------------
// -----------------------------Старая синхронизация без MPI---------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

//   // Для 3д случая добавлен цикл по z для синхронизации каждого яруса 2д границ \
//   а также отдельно без цикла добавлена синхронизация нижней и верхней граней \
//   относительно Oz.

//   // Для самого нижнего -1 и верхнего Nz + 1 яруса, необходимо сначала сам центр
//   // синхронизировать, а затем уже боковые узлы, так как они зависят от центра.

//   // ========================================================================

//   // ОБНОВЛЕНИЕ НИЖНЕЙ ГРАНИЦЫ ОТНОСИТЕЛЬНО Oz
//    // Отдельно угловые узлы в этом случае не обрабатываю, так как \
// если обновить бока (слева и справа), то для угловых узлов будет всё включено

//   // Центр нижней грани
// #pragma omp for collapse(2)
//   for (int64_t x = 0; x < Nx; ++x) {
//     for (int64_t y = 0; y < Ny; ++y) {
//       Ex(x, y, -1) = Ex(x, y, Nz - 1);
//       Ey(x, y, -1) = Ey(x, y, Nz - 1);
//       Ez(x, y, -1) = Ez(x, y, Nz - 1);
//       Bx(x, y, -1) = Bx(x, y, Nz - 1);
//       By(x, y, -1) = By(x, y, Nz - 1);
//       Bz(x, y, -1) = Bz(x, y, Nz - 1);
//     }
//   }
//   // Здесь вот не на 100% уверен насчёт Nz - 1, но вроде так
//   // Слева нижней грани (за исключением угловых узлов)
// #pragma omp for
//   for (int64_t y = 0; y < Ny; ++y) {
//     Ex(-1, y, -1) = Ex(Nx - 1, y, Nz - 1);
//     Ey(-1, y, -1) = Ey(Nx - 1, y, Nz - 1);
//     Ez(-1, y, -1) = Ez(Nx - 1, y, Nz - 1);
//     Bx(-1, y, -1) = Bx(Nx - 1, y, Nz - 1);
//     By(-1, y, -1) = By(Nx - 1, y, Nz - 1);
//     Bz(-1, y, -1) = Bz(Nx - 1, y, Nz - 1);
//   }
//   // Справа нижней грани (за исключением угловых узлов)
// #pragma omp for
//   for (int64_t y = 0; y < Ny; ++y) {
//     Ex(Nx, y, -1) = Ex(0, y, Nz - 1);
//     Ey(Nx, y, -1) = Ey(0, y, Nz - 1);
//     Ez(Nx, y, -1) = Ez(0, y, Nz - 1);
//     Bx(Nx, y, -1) = Bx(0, y, Nz - 1);
//     By(Nx, y, -1) = By(0, y, Nz - 1);
//     Bz(Nx, y, -1) = Bz(0, y, Nz - 1);
//   }
//   // Верх нижней грани (с учётом угловых узлов)
// #pragma omp for
//   for (int64_t x = -1; x < Nx + 1; ++x) {
//     Ex(x, Ny, -1) = Ex(x, 0, Nz - 1);
//     Ey(x, Ny, -1) = Ey(x, 0, Nz - 1);
//     Ez(x, Ny, -1) = Ez(x, 0, Nz - 1);
//     Bx(x, Ny, -1) = Bx(x, 0, Nz - 1);
//     By(x, Ny, -1) = By(x, 0, Nz - 1);
//     Bz(x, Ny, -1) = Bz(x, 0, Nz - 1);
//   }
//   // Низ нижней грани (с учётом угловых узлов)
// #pragma omp for
//   for (int64_t x = -1; x < Nx + 1; ++x) {
//     Ex(x, -1, -1) = Ex(x, Ny - 1, Nz - 1);
//     Ey(x, -1, -1) = Ey(x, Ny - 1, Nz - 1);
//     Ez(x, -1, -1) = Ez(x, Ny - 1, Nz - 1);
//     Bx(x, -1, -1) = Bx(x, Ny - 1, Nz - 1);
//     By(x, -1, -1) = By(x, Ny - 1, Nz - 1);
//     Bz(x, -1, -1) = Bz(x, Ny - 1, Nz - 1);
//   }
//   // ========================================================================
//   // ОБНОВЛЕНИЕ ВЕРХНЕЙ ГРАНИЦЫ ОТНОСИТЕЛЬНО Oz

//   // Центр верхней грани
// #pragma omp for collapse(2)
//   for (int64_t x = 0; x < Nx; ++x) {
//     for (int64_t y = 0; y < Ny; ++y) {
//       Ex(x, y, Nz) = Ex(x, y, 0);
//       Ey(x, y, Nz) = Ey(x, y, 0);
//       Ez(x, y, Nz) = Ez(x, y, 0);
//       Bx(x, y, Nz) = Bx(x, y, 0);
//       By(x, y, Nz) = By(x, y, 0);
//       Bz(x, y, Nz) = Bz(x, y, 0);
//     }
//   }
//   // Слева верхней грани (за исключением угловых узлов)
// #pragma omp for
//   for (int64_t y = 0; y < Ny; ++y) {
//     Ex(-1, y, Nz) = Ex(Nx - 1, y, 0);
//     Ey(-1, y, Nz) = Ey(Nx - 1, y, 0);
//     Ez(-1, y, Nz) = Ez(Nx - 1, y, 0);
//     Bx(-1, y, Nz) = Bx(Nx - 1, y, 0);
//     By(-1, y, Nz) = By(Nx - 1, y, 0);
//     Bz(-1, y, Nz) = Bz(Nx - 1, y, 0);
//   }
//   // Справа верхней грани (за исключением угловых узлов)
// #pragma omp for
//   for (int64_t y = 0; y < Ny; ++y) {
//     Ex(Nx, y, Nz) = Ex(0, y, 0);
//     Ey(Nx, y, Nz) = Ey(0, y, 0);
//     Ez(Nx, y, Nz) = Ez(0, y, 0);
//     Bx(Nx, y, Nz) = Bx(0, y, 0);
//     By(Nx, y, Nz) = By(0, y, 0);
//     Bz(Nx, y, Nz) = Bz(0, y, 0);
//   }
//   // Верх верхней грани (с учётом угловых узлов)
// #pragma omp for
//   for (int64_t x = -1; x < Nx + 1; ++x) {
//     Ex(x, Ny, Nz) = Ex(x, 0, 0);
//     Ey(x, Ny, Nz) = Ey(x, 0, 0);
//     Ez(x, Ny, Nz) = Ez(x, 0, 0);
//     Bx(x, Ny, Nz) = Bx(x, 0, 0);
//     By(x, Ny, Nz) = By(x, 0, 0);
//     Bz(x, Ny, Nz) = Bz(x, 0, 0);
//   }
//   // Низ верхней грани (с учётом угловых узлов)
// #pragma omp for
//   for (int64_t x = -1; x < Nx + 1; ++x) {
//     Ex(x, -1, Nz) = Ex(x, Ny - 1, 0);
//     Ey(x, -1, Nz) = Ey(x, Ny - 1, 0);
//     Ez(x, -1, Nz) = Ez(x, Ny - 1, 0);
//     Bx(x, -1, Nz) = Bx(x, Ny - 1, 0);
//     By(x, -1, Nz) = By(x, Ny - 1, 0);
//     Bz(x, -1, Nz) = Bz(x, Ny - 1, 0);
//   }

//   // ========================================================================
//   // Синхронизация границ для 2д случая
//   // for (int64_t z{-1}; z < Nz + 1; ++z) { странно. Вроде как должно быть 0, Nz
//   // Щас попробую так

//   // for (int64_t z{-1}; z < Nz + 1; ++z) {
// #pragma omp for
//   for (int64_t z = 0; z < Nz; ++z) {
//     // Слева 2D (за исключением угловых узлов)
//     for (int64_t y = 0; y < Ny; ++y) {
//       Ex(-1, y, z) = Ex(Nx - 1, y, z);
//       Ey(-1, y, z) = Ey(Nx - 1, y, z);
//       Ez(-1, y, z) = Ez(Nx - 1, y, z);
//       Bx(-1, y, z) = Bx(Nx - 1, y, z);
//       By(-1, y, z) = By(Nx - 1, y, z);
//       Bz(-1, y, z) = Bz(Nx - 1, y, z);
//     }
//     // Справа 2D (за исключением угловых узлов)
//     for (int64_t y = 0; y < Ny; ++y) {
//       Ex(Nx, y, z) = Ex(0, y, z);
//       Ey(Nx, y, z) = Ey(0, y, z);
//       Ez(Nx, y, z) = Ez(0, y, z);
//       Bx(Nx, y, z) = Bx(0, y, z);
//       By(Nx, y, z) = By(0, y, z);
//       Bz(Nx, y, z) = Bz(0, y, z);
//     }
//     // Верх 2D (с учётом угловых узлов)
//     for (int64_t x = -1; x < Nx + 1; ++x) {
//       Ex(x, Ny, z) = Ex(x, 0, z);
//       Ey(x, Ny, z) = Ey(x, 0, z);
//       Ez(x, Ny, z) = Ez(x, 0, z);
//       Bx(x, Ny, z) = Bx(x, 0, z);
//       By(x, Ny, z) = By(x, 0, z);
//       Bz(x, Ny, z) = Bz(x, 0, z);
//     }
//     // Низ 2D (с учётом угловых узлов)
//     for (int64_t x = -1; x < Nx + 1; ++x) {
//       Ex(x, -1, z) = Ex(x, Ny - 1, z);
//       Ey(x, -1, z) = Ey(x, Ny - 1, z);
//       Ez(x, -1, z) = Ez(x, Ny - 1, z);
//       Bx(x, -1, z) = Bx(x, Ny - 1, z);
//       By(x, -1, z) = By(x, Ny - 1, z);
//       Bz(x, -1, z) = Bz(x, Ny - 1, z);
//     }
//   }
}

std::ostream &FDTD::operator<<(std::ostream &os, const Component &comp) {
  switch (comp) {
  case Component::Ex:
    os << "Ex";
    break;
  case Component::Ey:
    os << "Ey";
    break;
  case Component::Ez:
    os << "Ez";
    break;
  case Component::Bx:
    os << "Bx";
    break;
  case Component::By:
    os << "By";
    break;
  case Component::Bz:
    os << "Bz";
    break;
  default:
    os << "Unknown component of the field (Ex, Ey, Ez, Bx, By, Bz)";
    break;
  }
  return os;
}
