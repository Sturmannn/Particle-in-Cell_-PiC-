#include "FDTD.hpp"
#include <string.h>

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
      dx{std::get<0>(dx_dy_dz)}, 
      dy{std::get<1>(dx_dy_dz)},
      dz{std::get<2>(dx_dy_dz)},
      dt{_dt} {
        // int world_size = 0;
        // MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        // if (world_size == 1) return;

        // // Определение общего размеры сетки для правильного определения dx, dy, dz

        // int rank = 0;
        // MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // int64_t global_Nx_Ny_Nz[3] = {0, 0, 0};
        // int64_t local_Nx_Ny_Nz[3] = {Nx, Ny, Nz};

        // if (rank == 0) {
        //   std::cout << "Nx = " << local_Nx_Ny_Nz[0] << " Ny = " << local_Nx_Ny_Nz[1] << " Nz = " << local_Nx_Ny_Nz[2] << '\n';
        // }

        // MPI_Allreduce(&local_Nx_Ny_Nz, &global_Nx_Ny_Nz, 3, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

        // if (rank == 0) {
        //   std::cout << "Nx = " << global_Nx_Ny_Nz[0] << " Ny = " << global_Nx_Ny_Nz[1] << " Nz = " << global_Nx_Ny_Nz[2] << '\n';
        // }

        // global_Nx_Ny_Nz[0] = std::max(global_Nx_Ny_Nz[0] - 2 * world_size, static_cast<int64_t>(1));
        // global_Nx_Ny_Nz[1] = std::max(global_Nx_Ny_Nz[1] - 2 * world_size, static_cast<int64_t>(1));
        // // Важно правильно настроить Nz, в случае, если сетка 2D
        // // ОПАСНОЕ МЕСТО. Не уверен, учтены ли здесь все случаи!!!
        // global_Nx_Ny_Nz[2] = std::max(global_Nx_Ny_Nz[2] - 2 * world_size, static_cast<int64_t>(1));

        // if (rank == 0) {
        //   std::cout << "Nx = " << global_Nx_Ny_Nz[0] << " Ny = " << global_Nx_Ny_Nz[1] << " Nz = " << global_Nx_Ny_Nz[2] << '\n';
        // }
        
        // dx = (get_bx() - get_ax()) / static_cast<double>(global_Nx_Ny_Nz[0]);
        // dy = (get_by() - get_ay()) / static_cast<double>(global_Nx_Ny_Nz[1]);
        // // Важно правильно настроить dz, в случае, если сетка 2D
        // dz = (global_Nx_Ny_Nz[2] > 1) ? (get_bz() - get_az()) / static_cast<double>(global_Nx_Ny_Nz[2]) : 0.0;

        // std::cout << "rank = " << rank << " dx = " << dx << " dy = " << dy << " dz = " << dz << '\n';

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
void FDTD::FDTD::field_update(const double t) {
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
    boundary_synchronization();
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
    boundary_synchronization();
  }
}
}

// Only for 2D
void FDTD::FDTD::field_update(const int64_t t) {
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
    boundary_synchronization();
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
    boundary_synchronization();
  }
}
}

void FDTD::FDTD::shifted_field_update(const double t) {}
void FDTD::FDTD::shifted_field_update(const int64_t t) {
  if (dt == 0.0) {
    std::cout << "Time step is null";
    exit(-1);
  }

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double B_dt = dt * 0.5;
  double E_dt = dt;

  int64_t _Nz = (Nz > 1) ? Nz : 0; // Учитываем случай 2D
  int64_t k_value = (_Nz > 1) ? 0 : -1; // Это эквивалентно -1 для 2D

  int64_t i{0};
  int64_t j{0};
  int64_t k{k_value};

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
      boundary_synchronization();
// #pragma omp for collapse(2)
      for (i = 0; i < Nx; ++i)
        for (j = 0; j < Ny; ++j) {
          Ex(i, j) = Ex(i, j) + C * E_dt * ((Bz(i, j) - Bz(i, j - 1)) / dy);
          Ey(i, j) = Ey(i, j) + C * E_dt * (-(Bz(i, j) - Bz(i - 1, j)) / dx);
          Ez(i, j) = Ez(i, j) + C * E_dt *
                                    ((By(i, j) - By(i - 1, j)) / dx -
                                     (Bx(i, j) - Bx(i, j - 1)) / dy);

        }
      boundary_synchronization();
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
      boundary_synchronization();
    }
  }
  }  /* else {
#pragma omp parallel private(i, j, k)
    {
      for (int64_t time = 0; time < t; ++time) // ПОПРАВИТЬ НА time += E_dt
      {
        // std::cout << "Hello from " << omp_get_thread_num() << std::endl;
#pragma omp for collapse(3)
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
        boundary_synchronization_3D();
#pragma omp for collapse(3)
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
        boundary_synchronization_3D();
#pragma omp for collapse(3)
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
        boundary_synchronization_3D();
      }
    }
  } */
}

// void FDTD::FDTD::shifted_field_update(const int64_t t) {
//   if (dt == 0.0) {
//     std::cout << "Time step is null";
//     exit(-1);
//   }
//   double B_dt = dt * 0.5;
//   double E_dt = dt;
//   int64_t i{0};
//   int64_t j{0};
//   int64_t k{0};
//   // double start = omp_get_wtime();

//   // #ifdef MPI

//   int rank = 0, mpi_comm_size = 0;
//   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   MPI_Comm_size(MPI_COMM_WORLD, &mpi_comm_size);
//   MPI_Status status;

//   int64_t Nx_local = get_Nx() + 2;
//   int64_t Ny_local{};
//   int64_t Nz_local = 1;

//   // Пока что хвост я запихиваю всё в root процесс (потом можно попробовать
//   // более равномерно распределить на другие процессы)
//   if (rank == 0) {
//     Ny_local = get_Ny() / mpi_comm_size + 2 + get_Ny() % mpi_comm_size;
//     std::cout << "Ny_local0 = " << Ny_local << std::endl;
//   } else {
//     Ny_local = get_Ny() / mpi_comm_size + 2;
//     std::cout << "Ny_local(drugie) = " << Ny_local << std::endl;
//   }

//   // if (mpi_comm_size == 1)
//   //{
//   //   Ny_local = (get_Ny() + 2) / mpi_comm_size;
//   //   std::cout << "Ny_local0 = " << Ny_local << std::endl;
//   // }
//   // else
//   //{
//   //   Ny_local = ((get_Ny() + 2) / mpi_comm_size) + 1; //   int64_t Ny_local
//   =
//   //   (get_Ny() + 2) / mpi_comm_size + 1; std::cout << "Ny_local(drugie) = "
//   <<
//   //   Ny_local << std::endl;
//   // }

//   // int64_t Nx_local = get_Nx();
//   // int64_t Ny_local = get_Ny() / mpi_comm_size; //   int64_t Ny_local =
//   // (get_Ny() + 2) / mpi_comm_size + 1; int64_t Nz_local = 1;

//   // Nx_local_size = get_Nx();
//   // Ny_local_size = get_Ny() / mpi_comm_size;

//   // if (rank == 0)
//   //{
//   //   for (size_t i = 0; i < 144; ++i)
//   //   {
//   //     *(Bx.data() + i) = i;
//   //     //std::cout << "StartBx[" << i << "] = " << *(Bx.data() + i) <<
//   //     std::endl;
//   //   }
//   // }

//   Field::ComputingField Ex_local(Nx_local - 2, Ny_local - 2);
//   Field::ComputingField Ey_local(Nx_local - 2, Ny_local - 2);
//   Field::ComputingField Ez_local(Nx_local - 2, Ny_local - 2);
//   Field::ComputingField Bx_local(Nx_local - 2, Ny_local - 2);
//   Field::ComputingField By_local(Nx_local - 2, Ny_local - 2);
//   Field::ComputingField Bz_local(Nx_local - 2, Ny_local - 2);
//   // std::cout << "Ex_local[" << rank << "] = " << Nx_local - 2 << " " <<
//   // Ny_local - 2 << std::endl;
//   std::cout << "Ex_local[" << rank << "] = " << Ex_local.fullsize()
//             << std::endl;
//   std::cout << "Ex[" << rank << "] = " << Ex.fullsize() << std::endl;

//   std::vector<int> sendcounts{};
//   int recvcountGather = Nx_local * Ny_local;
//   if (rank == 0)
//     sendcounts.resize(mpi_comm_size);
//   // Собираю на root процессе размер данных для рассылки
//   MPI_Gather(&recvcountGather, 1, MPI_INT, sendcounts.data(), 1, MPI_INT, 0,
//              MPI_COMM_WORLD);

//   if (rank == 0) {
//     for (int i = 0; i < mpi_comm_size; ++i)
//       std::cout << "sendconts[" << i << "] = " << sendcounts[i] << std::endl;
//   }
//   std::vector<int> displs(mpi_comm_size);
//   // Смещения считаются только на отсылающем процессе
//   if (rank == 0) {
//     for (int i = 0; i < displs.size(); ++i) {
//       displs[i] =
//           (static_cast<int>(Ny_local) - 2) * static_cast<int>(Nx_local) * i;
//       std::cout << "displs[" << i << "] = " << displs[i] << std::endl;
//     }
//   }

//   MPI_Scatterv(Ex.data(), sendcounts.data(), displs.data(), MPI_DOUBLE,
//                Ex_local.data(), static_cast<int>(Ex_local.fullsize()),
//                MPI_DOUBLE, 0, MPI_COMM_WORLD);
//   MPI_Scatterv(Ey.data(), sendcounts.data(), displs.data(), MPI_DOUBLE,
//                Ey_local.data(), static_cast<int>(Ey_local.fullsize()),
//                MPI_DOUBLE, 0, MPI_COMM_WORLD);
//   MPI_Scatterv(Ez.data(), sendcounts.data(), displs.data(), MPI_DOUBLE,
//                Ez_local.data(), static_cast<int>(Ez_local.fullsize()),
//                MPI_DOUBLE, 0, MPI_COMM_WORLD);

//   MPI_Scatterv(Bx.data(), sendcounts.data(), displs.data(), MPI_DOUBLE,
//                Bx_local.data(), static_cast<int>(Bx_local.fullsize()),
//                MPI_DOUBLE, 0, MPI_COMM_WORLD);
//   MPI_Scatterv(By.data(), sendcounts.data(), displs.data(), MPI_DOUBLE,
//                By_local.data(), static_cast<int>(By_local.fullsize()),
//                MPI_DOUBLE, 0, MPI_COMM_WORLD);
//   MPI_Scatterv(Bz.data(), sendcounts.data(), displs.data(), MPI_DOUBLE,
//                Bz_local.data(), static_cast<int>(Bz_local.fullsize()),
//                MPI_DOUBLE, 0, MPI_COMM_WORLD);

//   // if (rank == 0)
//   // for (int i = 0; i < Bx.size(); ++i)
//   //{
//   //   *(Bx.data() + i) = 0;
//   // }
//   //
//   // if (rank == 0)
//   //{
//   //   for (int i = 0; i < By_local.size(); ++i)
//   //   {
//   //     std::cout << "0By[" << i << "] = " << *(By_local.data() + i) <<
//   //     std::endl;
//   //   }
//   // }

//   for (int64_t time = 0ull; time < t; time++) {
//     MPI_Barrier(MPI_COMM_WORLD);
//     // for (i = 0; i < Nx_local; ++i)
//     //   for (j = 0; j < Ny_local; ++j)
//     //     for (k = 0; k < Nz_local; ++k)
//     for (i = 0; i < Nx_local - 2; ++i)
//       for (j = 0; j < Ny_local - 2; ++j)
//         for (k = 0; k < Nz_local; ++k) {
//           // std::cout << Ez_local(i + 1, j, k) << " " << Ez_local(i, j, k)
//           << "
//           // " << (Ez_local(i + 1, j, k) - Ez_local(i, j, k)) << std::endl;
//           // if (By_local(i, j, k) != 0)
//           //{
//           //   std::cout << i << " " << j << " " << k << " " << time <<
//           //   std::endl;
//           // }
//           Bx_local(i, j, k) =
//               Bx_local(i, j, k) +
//               C * B_dt *
//                   (/*(Ey_local(i, j, k + 1) - Ey_local(i, j, k)) / dz*/ -(
//                        Ez_local(i, j + 1, k) - Ez_local(i, j, k)) /
//                    dy);
//           By_local(i, j, k) =
//               By_local(i, j, k) +
//               C * B_dt *
//                   ((Ez_local(i + 1, j, k) - Ez_local(i, j, k)) /
//                    dx /*- (Ex_local(i, j, k + 1) - Ex_local(i, j, k)) /
//                    dz*/);
//           Bz_local(i, j, k) =
//               Bz_local(i, j, k) +
//               C * B_dt *
//                   ((Ex_local(i, j + 1, k) - Ex_local(i, j, k)) / dy -
//                    (Ey_local(i + 1, j, k) - Ey_local(i, j, k)) / dx);
//         }

//     // Обновляю бока
//     for (int64_t y = 0; y < Bx_local.get_Ny(); ++y) {
//       Bx_local(-1, y) = Bx_local(Bx_local.get_Nx() - 1, y); // слева
//       // Bx_local(Bx_local.get_Nx(), y) = Bx_local(0, y); // справа

//       By_local(-1, y) = By_local(By_local.get_Nx() - 1, y); // слева
//       By_local(By_local.get_Nx(), y) = By_local(0, y);      // справа

//       Bz_local(-1, y) = Bz_local(Bz_local.get_Nx() - 1, y); // слева
//       Bz_local(Bz_local.get_Nx(), y) = Bz_local(0, y);      // справа
//     }

//     // Обновляю границу - верх 2-х мерной системы
//     if (rank == 0 && mpi_comm_size != 1) {
//       // Получаю верхнюю строчку последнего процесса для синхронизации
//       нулевой
//       // строчки нулевого процесса Замечание: здесь get_Nx() = Nx_local - 2
//       MPI_Recv(Bx_local.data() + 1, static_cast<int>(Bx_local.get_Nx()),
//                MPI_DOUBLE, mpi_comm_size - 1, 0, MPI_COMM_WORLD,
//                &status); // снизу
//       MPI_Recv(By_local.data() + 1, static_cast<int>(By_local.get_Nx()),
//                MPI_DOUBLE, mpi_comm_size - 1, 1, MPI_COMM_WORLD,
//                &status); // снизу
//       MPI_Recv(Bz_local.data() + 1, static_cast<int>(Bz_local.get_Nx()),
//                MPI_DOUBLE, mpi_comm_size - 1, 2, MPI_COMM_WORLD,
//                &status); // снизу

//       MPI_Recv(Bx_local.data(), 1, MPI_DOUBLE, mpi_comm_size - 1, 3,
//                MPI_COMM_WORLD, &status); // левый нижний узел
//       MPI_Recv(By_local.data(), 1, MPI_DOUBLE, mpi_comm_size - 1, 4,
//                MPI_COMM_WORLD, &status); // левый нижний узел
//       MPI_Recv(Bz_local.data(), 1, MPI_DOUBLE, mpi_comm_size - 1, 5,
//                MPI_COMM_WORLD, &status); // левый нижний узел

//       MPI_Recv(Bx_local.data() + static_cast<int>(Bx_local.get_Nx()) + 1, 1,
//                MPI_DOUBLE, mpi_comm_size - 1, 6, MPI_COMM_WORLD,
//                &status); // правый нижний узел
//       MPI_Recv(By_local.data() + static_cast<int>(By_local.get_Nx()) + 1, 1,
//                MPI_DOUBLE, mpi_comm_size - 1, 7, MPI_COMM_WORLD,
//                &status); // правый нижний узел
//       MPI_Recv(Bz_local.data() + static_cast<int>(Bz_local.get_Nx()) + 1, 1,
//                MPI_DOUBLE, mpi_comm_size - 1, 8, MPI_COMM_WORLD,
//                &status); // правый нижний узел

//       MPI_Ssend(&Bx_local(0, 0), static_cast<int>(Bx_local.get_Nx()),
//                 MPI_DOUBLE, mpi_comm_size - 1, 9,
//                 MPI_COMM_WORLD); // Верх последнего процесса
//       MPI_Ssend(&By_local(0, 0), static_cast<int>(By_local.get_Nx()),
//                 MPI_DOUBLE, mpi_comm_size - 1, 10,
//                 MPI_COMM_WORLD); // Верх последнего процесса
//       MPI_Ssend(&Bz_local(0, 0), static_cast<int>(Bz_local.get_Nx()),
//                 MPI_DOUBLE, mpi_comm_size - 1, 11,
//                 MPI_COMM_WORLD); // Верх последнего процесса

//       MPI_Ssend(&Bx_local(0, 0), 1, MPI_DOUBLE, mpi_comm_size - 1, 12,
//                 MPI_COMM_WORLD); // левый нижний угол (для правого врехнего
//                 угла
//                                  // посл. процесса)
//       MPI_Ssend(&By_local(0, 0), 1, MPI_DOUBLE, mpi_comm_size - 1, 13,
//                 MPI_COMM_WORLD); // левый нижний угол (для правого врехнего
//                 угла
//                                  // посл. процесса)
//       MPI_Ssend(&Bz_local(0, 0), 1, MPI_DOUBLE, mpi_comm_size - 1, 14,
//                 MPI_COMM_WORLD); // левый нижний угол (для правого врехнего
//                 угла
//                                  // посл. процесса)

//       MPI_Ssend(&Bx_local(Bx_local.get_Nx() - 1, 0), 1, MPI_DOUBLE,
//                 mpi_comm_size - 1, 15,
//                 MPI_COMM_WORLD); // правый нижний угол (для лев. верх. угла
//                                  // посл. процесса)
//       MPI_Ssend(&By_local(By_local.get_Nx() - 1, 0), 1, MPI_DOUBLE,
//                 mpi_comm_size - 1, 16,
//                 MPI_COMM_WORLD); // правый нижний угол (для лев. верх. угла
//                                  // посл. процесса)
//       MPI_Ssend(&Bz_local(Bz_local.get_Nx() - 1, 0), 1, MPI_DOUBLE,
//                 mpi_comm_size - 1, 17,
//                 MPI_COMM_WORLD); // правый нижний угол (для лев. верх. угла
//                                  // посл. процесса)
//     }
//     if (rank == mpi_comm_size - 1 && mpi_comm_size != 1) {
//       // for (int64_t y = 0; y < Bx_local.get_Ny(); ++y)
//       //{
//       //   Bx_local(-1, y) = Bx_local(Bx_local.get_Nx() - 1, y); // слева
//       //   Bx_local(Bx_local.get_Nx(), y) = Bx_local(0, y); // справа
//       //
//       //   By_local(-1, y) = By_local(By_local.get_Nx() - 1, y); // слева
//       //   By_local(By_local.get_Nx(), y) = By_local(0, y); // справа
//       //
//       //   Bz_local(-1, y) = Bz_local(Bz_local.get_Nx() - 1, y); // слева
//       //   Bz_local(Bz_local.get_Nx(), y) = Bz_local(0, y); // справа
//       // }
//       //  Отправляю верхнюю строчку последнего процесса нулевому процессу для
//       //  нулевой строчки

//       MPI_Ssend(Bx_local.data() + Nx_local * (Ny_local - 2) + 1,
//                 static_cast<int>(Bx_local.get_Nx()), MPI_DOUBLE, 0, 0,
//                 MPI_COMM_WORLD); // Низ нулевого процесса (тэг 0)
//       MPI_Ssend(By_local.data() + Nx_local * (Ny_local - 2) + 1,
//                 static_cast<int>(By_local.get_Nx()), MPI_DOUBLE, 0, 1,
//                 MPI_COMM_WORLD); // Низ нулевого процесса (тэг 1)
//       MPI_Ssend(Bz_local.data() + Nx_local * (Ny_local - 2) + 1,
//                 static_cast<int>(Bz_local.get_Nx()), MPI_DOUBLE, 0, 2,
//                 MPI_COMM_WORLD); // Низ нулевого процесса (тэг 2)

//       MPI_Ssend(Bx_local.data() + Nx_local * (Ny_local - 1) - 2, 1,
//       MPI_DOUBLE,
//                 0, 3,
//                 MPI_COMM_WORLD); // левый нижний узел нулевого процесса (тэг
//                 3)
//       MPI_Ssend(By_local.data() + Nx_local * (Ny_local - 1) - 2, 1,
//       MPI_DOUBLE,
//                 0, 4,
//                 MPI_COMM_WORLD); // левый нижний узел нулевого процесса (тэг
//                 4)
//       MPI_Ssend(Bz_local.data() + Nx_local * (Ny_local - 1) - 2, 1,
//       MPI_DOUBLE,
//                 0, 5,
//                 MPI_COMM_WORLD); // левый нижний узел нулевого процесса (тэг
//                 5)

//       MPI_Ssend(Bx_local.data() + Nx_local * (Ny_local - 2) + 1, 1,
//       MPI_DOUBLE,
//                 0, 6,
//                 MPI_COMM_WORLD); // правый нижний узел узел нулевого процесса
//       MPI_Ssend(By_local.data() + Nx_local * (Ny_local - 2) + 1, 1,
//       MPI_DOUBLE,
//                 0, 7,
//                 MPI_COMM_WORLD); // правый нижний узел узел нулевого процесса
//       MPI_Ssend(Bz_local.data() + Nx_local * (Ny_local - 2) + 1, 1,
//       MPI_DOUBLE,
//                 0, 8,
//                 MPI_COMM_WORLD); // правый нижний узел узел нулевого процесса

//       MPI_Recv(&Bx_local(0, Bx_local.get_Ny()),
//                static_cast<int>(Bx_local.get_Nx()), MPI_DOUBLE, 0, 9,
//                MPI_COMM_WORLD,
//                &status); // сверху (получение нулевой строки нулевого
//                процесса)
//       MPI_Recv(&By_local(0, By_local.get_Ny()),
//                static_cast<int>(By_local.get_Nx()), MPI_DOUBLE, 0, 10,
//                MPI_COMM_WORLD,
//                &status); // сверху (получение нулевой строки нулевого
//                процесса)
//       MPI_Recv(&Bz_local(0, Bz_local.get_Ny()),
//                static_cast<int>(Bz_local.get_Nx()), MPI_DOUBLE, 0, 11,
//                MPI_COMM_WORLD,
//                &status); // сверху (получение нулевой строки нулевого
//                процесса)

//       MPI_Recv(&Bx_local(Bx_local.get_Nx(), Bx_local.get_Ny()), 1,
//       MPI_DOUBLE,
//                0, 12, MPI_COMM_WORLD,
//                &status); // правый верхний узел (получение из нулевого
//                процесса)
//       MPI_Recv(&By_local(By_local.get_Nx(), By_local.get_Ny()), 1,
//       MPI_DOUBLE,
//                0, 13, MPI_COMM_WORLD,
//                &status); // правый верхний узел (получение из нулевого
//                процесса)
//       MPI_Recv(&Bz_local(Bz_local.get_Nx(), Bz_local.get_Ny()), 1,
//       MPI_DOUBLE,
//                0, 14, MPI_COMM_WORLD,
//                &status); // правый верхний узел (получение из нулевого
//                процесса)

//       MPI_Recv(&Bx_local(-1, Bx_local.get_Ny()), 1, MPI_DOUBLE, 0, 15,
//                MPI_COMM_WORLD,
//                &status); // левый верхний узел (получение из нулевого
//                процесса)
//       MPI_Recv(&By_local(-1, By_local.get_Ny()), 1, MPI_DOUBLE, 0, 16,
//                MPI_COMM_WORLD,
//                &status); // левый верхний узел (получение из нулевого
//                процесса)
//       MPI_Recv(&Bz_local(-1, Bz_local.get_Ny()), 1, MPI_DOUBLE, 0, 17,
//                MPI_COMM_WORLD,
//                &status); // левый верхний узел (получение из нулевого
//                процесса)
//     }
//     MPI_Barrier(MPI_COMM_WORLD);
//     // for (size_t i = 0; i < 84; ++i)
//     //   std::cout << "11Bx[" << i << "] = " << *(Bx_local.data() + i) <<
//     "\n";

//     if (mpi_comm_size == 1) {
//       for (int64_t x = 0; x < Bx_local.get_Nx(); ++x) {
//         *(Bx_local.data() + x + 1) =
//             Bx_local(x, Bx_local.get_Ny() - 1);          // снизу
//         Bx_local(x, Bx_local.get_Ny()) = Bx_local(x, 0); // сверху

//         *(By_local.data() + x + 1) =
//             By_local(x, By_local.get_Ny() - 1);          // снизу
//         By_local(x, By_local.get_Ny()) = By_local(x, 0); // сверху

//         *(Bz_local.data() + x + 1) =
//             Bz_local(x, Bz_local.get_Ny() - 1);          // снизу
//         Bz_local(x, Bz_local.get_Ny()) = Bz_local(x, 0); // сверху
//       }
//       for (int64_t y = 0; y < Bx_local.get_Ny(); ++y) {
//         Bx_local(-1, y) = Bx_local(Bx_local.get_Nx() - 1, y); // слева
//         // Bx_local(Bx_local.get_Nx(), y) = Bx_local(0, y); // справа

//         By_local(-1, y) = By_local(By_local.get_Nx() - 1, y); // слева
//         By_local(By_local.get_Nx(), y) = By_local(0, y);      // справа

//         Bz_local(-1, y) = Bz_local(Bz_local.get_Nx() - 1, y); // слева
//         Bz_local(Bz_local.get_Nx(), y) = Bz_local(0, y);      // справа
//       }
//       Bx_local(-1, -1) = Bx_local(Bx_local.get_Nx() - 1,
//                                   Bx_local.get_Ny() - 1); // левый нижний
//                                   узел
//       Bx_local(-1, Bx_local.get_Ny()) =
//           Bx_local(Bx_local.get_Nx() - 1, 0); // левый верхний узел
//       Bx_local(Bx_local.get_Nx(), Bx_local.get_Ny()) =
//           Bx_local(0, 0); // правый верхний узел
//       *(Bx_local.data() + Bx_local.get_Nx() + 1) =
//           Bx_local(0, Bx_local.get_Ny() - 1); // правый нижний узел

//       By_local(-1, -1) = By_local(By_local.get_Nx() - 1,
//                                   By_local.get_Ny() - 1); // левый нижний
//                                   узел
//       By_local(-1, By_local.get_Ny()) =
//           By_local(By_local.get_Nx() - 1, 0); // левый верхний узел
//       By_local(By_local.get_Nx(), By_local.get_Ny()) =
//           By_local(0, 0); // правый верхний узел
//       *(By_local.data() + By_local.get_Nx() + 1) =
//           By_local(0, By_local.get_Ny() - 1); // правый нижний узел

//       Bz_local(-1, -1) = Bz_local(Bz_local.get_Nx() - 1,
//                                   Bz_local.get_Ny() - 1); // левый нижний
//                                   узел
//       Bz_local(-1, Bz_local.get_Ny()) =
//           Bz_local(Bz_local.get_Nx() - 1, 0); // левый верхний узел
//       Bz_local(Bz_local.get_Nx(), Bz_local.get_Ny()) =
//           Bz_local(0, 0); // правый верхний узел
//       *(Bz_local.data() + Bz_local.get_Nx() + 1) =
//           Bz_local(0, Bz_local.get_Ny() - 1); // правый нижний узел
//     }

//     // if (rank == 0)
//     //{
//     //       for (size_t i = 0; i < Nx_local; ++i)
//     //       {
//     //         *(Bx_local.data() + 60 + i) = 300 + i;
//     //       }
//     // }
//     // if (rank == 1)
//     //{
//     //     for (size_t i = 0; i < Nx_local; ++i)
//     //     {
//     //       *(Bx_local.data() + 12 + i) = 400 + i;
//     //     }
//     //
//     // }

//     // Синхронизация смежных строк (для 2-мерного случая)
//     if (mpi_comm_size > 1) {
//       if (rank != mpi_comm_size - 1 && rank != 0) {
//         MPI_Sendrecv(Bx_local.data() + Nx_local * (Ny_local - 2),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank,
//                      Bx_local.data() + Nx_local * (Ny_local - 1),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank +
//                      1, MPI_COMM_WORLD, &status); // сверху получил-отправил
//         MPI_Sendrecv(By_local.data() + Nx_local * (Ny_local - 2),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank,
//                      By_local.data() + Nx_local * (Ny_local - 1),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank +
//                      1, MPI_COMM_WORLD, &status); // сверху получил-отправил
//         MPI_Sendrecv(Bz_local.data() + Nx_local * (Ny_local - 2),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank,
//                      Bz_local.data() + Nx_local * (Ny_local - 1),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank +
//                      1, MPI_COMM_WORLD, &status); // сверху получил-отправил

//         MPI_Sendrecv(Bx_local.data() + Nx_local, static_cast<int>(Nx_local),
//                      MPI_DOUBLE, rank - 1, rank, Bx_local.data(),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank - 1, rank -
//                      1, MPI_COMM_WORLD, &status); // снизу получил-отправил
//         MPI_Sendrecv(By_local.data() + Nx_local, static_cast<int>(Nx_local),
//                      MPI_DOUBLE, rank - 1, rank, By_local.data(),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank - 1, rank -
//                      1, MPI_COMM_WORLD, &status); // снизу получил-отправил
//         MPI_Sendrecv(Bz_local.data() + Nx_local, static_cast<int>(Nx_local),
//                      MPI_DOUBLE, rank - 1, rank, Bz_local.data(),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank - 1, rank -
//                      1, MPI_COMM_WORLD, &status); // снизу получил-отправил
//       } else if (rank == 0) {
//         MPI_Sendrecv(Bx_local.data() + Nx_local * (Ny_local - 2),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank,
//                      Bx_local.data() + Nx_local * (Ny_local - 1),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank +
//                      1, MPI_COMM_WORLD, &status); // сверху получил-отправил
//         MPI_Sendrecv(By_local.data() + Nx_local * (Ny_local - 2),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank,
//                      By_local.data() + Nx_local * (Ny_local - 1),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank +
//                      1, MPI_COMM_WORLD, &status); // сверху получил-отправил
//         MPI_Sendrecv(Bz_local.data() + Nx_local * (Ny_local - 2),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank,
//                      Bz_local.data() + Nx_local * (Ny_local - 1),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank +
//                      1, MPI_COMM_WORLD, &status); // сверху получил-отправил
//       } else                                   // if (rank == mpi_comm_size -
//       1)
//       {

//         MPI_Sendrecv(Bx_local.data() + Nx_local, static_cast<int>(Nx_local),
//                      MPI_DOUBLE, rank - 1, rank, Bx_local.data(),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank - 1, rank -
//                      1, MPI_COMM_WORLD, &status); // снизу получил-отправил
//         MPI_Sendrecv(By_local.data() + Nx_local, static_cast<int>(Nx_local),
//                      MPI_DOUBLE, rank - 1, rank, By_local.data(),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank - 1, rank -
//                      1, MPI_COMM_WORLD, &status); // снизу получил-отправил
//         MPI_Sendrecv(Bz_local.data() + Nx_local, static_cast<int>(Nx_local),
//                      MPI_DOUBLE, rank - 1, rank, Bz_local.data(),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank - 1, rank -
//                      1, MPI_COMM_WORLD, &status); // снизу получил-отправил
//       }
//       MPI_Barrier(MPI_COMM_WORLD);
//       // if (rank == 0)
//       //{
//       //   for (size_t i = 0; i < 84; ++i)
//       //     std::cout << "0Bx[" << i << "] = " << *(Bx_local.data() + i) <<
//       //     "\n";
//       // }
//       // if (rank == 1)
//       //{
//       //   for (size_t i = 0; i < 84; ++i)
//       //     std::cout << "1Bx[" << i << "] = " << *(Bx_local.data() + i) <<
//       //     "\n";
//       // }
//     }
//     //
//     ==============================================================================================================================================

//     // for (i = 0ull; i < Nx_local; ++i)
//     //   for (j = 0ull; j < Ny_local; ++j)
//     //     for (k = 0ull; k < Nz_local; ++k)
//     MPI_Barrier(MPI_COMM_WORLD);
//     for (i = 0; i < Nx_local - 2; ++i)
//       for (j = 0; j < Ny_local - 2; ++j)
//         for (k = 0; k < Nz_local; ++k) {
//           Ex_local(i, j, k) =
//               Ex_local(i, j, k) +
//               C * E_dt *
//                   ((Bz_local(i, j, k) - Bz_local(i, j - 1, k)) /
//                    dy /*- (By_local(i, j, k) - By_local(i, j, k - 1)) /
//                    dz*/);
//           Ey_local(i, j, k) =
//               Ey_local(i, j, k) +
//               C * E_dt *
//                   (/*(Bx_local(i, j, k) - Bx_local(i, j, k - 1)) / dz*/ -(
//                        Bz_local(i, j, k) - Bz_local(i - 1, j, k)) /
//                    dx);
//           Ez_local(i, j, k) =
//               Ez_local(i, j, k) +
//               C * E_dt *
//                   ((By_local(i, j, k) - By_local(i - 1, j, k)) / dx -
//                    (Bx_local(i, j, k) - Bx_local(i, j - 1, k)) / dy);
//         }
//     // Обновляю бока
//     for (int64_t y = 0; y < Ex_local.get_Ny(); ++y) {
//       Ex_local(-1, y) = Ex_local(Ex_local.get_Nx() - 1, y); // слева
//       Ex_local(Ex_local.get_Nx(), y) = Ex_local(0, y);      // справа

//       Ey_local(-1, y) = Ey_local(Ey_local.get_Nx() - 1, y); // слева
//       Ey_local(Ey_local.get_Nx(), y) = Ey_local(0, y);      // справа

//       Ez_local(-1, y) = Ez_local(Ez_local.get_Nx() - 1, y); // слева
//       Ez_local(Ez_local.get_Nx(), y) = Ez_local(0, y);      // справа
//     }

//     // Обновляю границу - верх 2-х мерной системы
//     if (rank == 0 && mpi_comm_size != 1) {

//       // for (int64_t y = 0; y < Ex_local.get_Ny(); ++y)
//       //{
//       //   Ex_local(-1, y) = Ex_local(Ex_local.get_Nx() - 1, y); // слева
//       //   Ex_local(Ex_local.get_Nx(), y) = Ex_local(0, y); // справа
//       //
//       //   Ey_local(-1, y) = Ey_local(Ey_local.get_Nx() - 1, y); // слева
//       //   Ey_local(Ey_local.get_Nx(), y) = Ey_local(0, y); // справа
//       //
//       //   Ez_local(-1, y) = Ez_local(Ez_local.get_Nx() - 1, y); // слева
//       //   Ez_local(Ez_local.get_Nx(), y) = Ez_local(0, y); // справа
//       // }
//       //  Получаю верхнюю строчку последнего процесса для синхронизации
//       нулевой
//       //  строчки нулевого процесса Замечание: здесь get_Nx() = Nx_local - 2
//       MPI_Recv(Ex_local.data() + 1, static_cast<int>(Ex_local.get_Nx()),
//                MPI_DOUBLE, mpi_comm_size - 1, 0, MPI_COMM_WORLD,
//                &status); // снизу
//       MPI_Recv(Ey_local.data() + 1, static_cast<int>(Ey_local.get_Nx()),
//                MPI_DOUBLE, mpi_comm_size - 1, 1, MPI_COMM_WORLD,
//                &status); // снизу
//       MPI_Recv(Ez_local.data() + 1, static_cast<int>(Ez_local.get_Nx()),
//                MPI_DOUBLE, mpi_comm_size - 1, 2, MPI_COMM_WORLD,
//                &status); // снизу

//       MPI_Recv(Ex_local.data(), 1, MPI_DOUBLE, mpi_comm_size - 1, 3,
//                MPI_COMM_WORLD, &status); // левый нижний узел
//       MPI_Recv(Ey_local.data(), 1, MPI_DOUBLE, mpi_comm_size - 1, 4,
//                MPI_COMM_WORLD, &status); // левый нижний узел
//       MPI_Recv(Ez_local.data(), 1, MPI_DOUBLE, mpi_comm_size - 1, 5,
//                MPI_COMM_WORLD, &status); // левый нижний узел

//       MPI_Recv(Ex_local.data() + static_cast<int>(Ex_local.get_Nx()) + 1, 1,
//                MPI_DOUBLE, mpi_comm_size - 1, 6, MPI_COMM_WORLD,
//                &status); // правый нижний узел
//       MPI_Recv(Ey_local.data() + static_cast<int>(Ey_local.get_Nx()) + 1, 1,
//                MPI_DOUBLE, mpi_comm_size - 1, 7, MPI_COMM_WORLD,
//                &status); // правый нижний узел
//       MPI_Recv(Ez_local.data() + static_cast<int>(Ez_local.get_Nx()) + 1, 1,
//                MPI_DOUBLE, mpi_comm_size - 1, 8, MPI_COMM_WORLD,
//                &status); // правый нижний узел

//       MPI_Ssend(&Ex_local(0, 0), static_cast<int>(Ex_local.get_Nx()),
//                 MPI_DOUBLE, mpi_comm_size - 1, 9,
//                 MPI_COMM_WORLD); // Верх последнего процесса
//       MPI_Ssend(&Ey_local(0, 0), static_cast<int>(Ey_local.get_Nx()),
//                 MPI_DOUBLE, mpi_comm_size - 1, 10,
//                 MPI_COMM_WORLD); // Верх последнего процесса
//       MPI_Ssend(&Ez_local(0, 0), static_cast<int>(Ez_local.get_Nx()),
//                 MPI_DOUBLE, mpi_comm_size - 1, 11,
//                 MPI_COMM_WORLD); // Верх последнего процесса

//       MPI_Ssend(&Ex_local(0, 0), 1, MPI_DOUBLE, mpi_comm_size - 1, 12,
//                 MPI_COMM_WORLD); // левый нижний угол (для правого врехнего
//                 угла
//                                  // посл. процесса)
//       MPI_Ssend(&Ey_local(0, 0), 1, MPI_DOUBLE, mpi_comm_size - 1, 13,
//                 MPI_COMM_WORLD); // левый нижний угол (для правого врехнего
//                 угла
//                                  // посл. процесса)
//       MPI_Ssend(&Ez_local(0, 0), 1, MPI_DOUBLE, mpi_comm_size - 1, 14,
//                 MPI_COMM_WORLD); // левый нижний угол (для правого врехнего
//                 угла
//                                  // посл. процесса)

//       MPI_Ssend(&Ex_local(Ex_local.get_Nx() - 1, 0), 1, MPI_DOUBLE,
//                 mpi_comm_size - 1, 15,
//                 MPI_COMM_WORLD); // правый нижний угол (для лев. верх. угла
//                                  // посл. процесса)
//       MPI_Ssend(&Ey_local(Ey_local.get_Nx() - 1, 0), 1, MPI_DOUBLE,
//                 mpi_comm_size - 1, 16,
//                 MPI_COMM_WORLD); // правый нижний угол (для лев. верх. угла
//                                  // посл. процесса)
//       MPI_Ssend(&Ez_local(Ez_local.get_Nx() - 1, 0), 1, MPI_DOUBLE,
//                 mpi_comm_size - 1, 17,
//                 MPI_COMM_WORLD); // правый нижний угол (для лев. верх. угла
//                                  // посл. процесса)
//     }
//     if (rank == mpi_comm_size - 1 && mpi_comm_size != 1) {

//       // for (int64_t y = 0; y < Ex_local.get_Ny(); ++y)
//       //{
//       //   Ex_local(-1, y) = Ex_local(Ex_local.get_Nx() - 1, y); // слева
//       //   Ex_local(Ex_local.get_Nx(), y) = Ex_local(0, y); // справа
//       //
//       //   Ey_local(-1, y) = Ey_local(Ey_local.get_Nx() - 1, y); // слева
//       //   Ey_local(Ey_local.get_Nx(), y) = Ey_local(0, y); // справа
//       //
//       //   Ez_local(-1, y) = Ez_local(Ez_local.get_Nx() - 1, y); // слева
//       //   Ez_local(Ez_local.get_Nx(), y) = Ez_local(0, y); // справа
//       // }

//       // Отправляю верхнюю строчку последнего процесса нулевому процессу для
//       // нулевой строчки

//       MPI_Ssend(Ex_local.data() + Nx_local * (Ny_local - 2) + 1,
//                 static_cast<int>(Ex_local.get_Nx()), MPI_DOUBLE, 0, 0,
//                 MPI_COMM_WORLD); // Низ нулевого процесса (тэг 0)
//       MPI_Ssend(Ey_local.data() + Nx_local * (Ny_local - 2) + 1,
//                 static_cast<int>(Ey_local.get_Nx()), MPI_DOUBLE, 0, 1,
//                 MPI_COMM_WORLD); // Низ нулевого процесса (тэг 1)
//       MPI_Ssend(Ez_local.data() + Nx_local * (Ny_local - 2) + 1,
//                 static_cast<int>(Ez_local.get_Nx()), MPI_DOUBLE, 0, 2,
//                 MPI_COMM_WORLD); // Низ нулевого процесса (тэг 2)

//       MPI_Ssend(Ex_local.data() + Nx_local * (Ny_local - 1) - 2, 1,
//       MPI_DOUBLE,
//                 0, 3,
//                 MPI_COMM_WORLD); // левый нижний узел нулевого процесса (тэг
//                 3)
//       MPI_Ssend(Ey_local.data() + Nx_local * (Ny_local - 1) - 2, 1,
//       MPI_DOUBLE,
//                 0, 4,
//                 MPI_COMM_WORLD); // левый нижний узел нулевого процесса (тэг
//                 4)
//       MPI_Ssend(Ez_local.data() + Nx_local * (Ny_local - 1) - 2, 1,
//       MPI_DOUBLE,
//                 0, 5,
//                 MPI_COMM_WORLD); // левый нижний узел нулевого процесса (тэг
//                 5)

//       MPI_Ssend(Ex_local.data() + Nx_local * (Ny_local - 2) + 1, 1,
//       MPI_DOUBLE,
//                 0, 6,
//                 MPI_COMM_WORLD); // правый нижний узел узел нулевого процесса
//       MPI_Ssend(Ey_local.data() + Nx_local * (Ny_local - 2) + 1, 1,
//       MPI_DOUBLE,
//                 0, 7,
//                 MPI_COMM_WORLD); // правый нижний узел узел нулевого процесса
//       MPI_Ssend(Ez_local.data() + Nx_local * (Ny_local - 2) + 1, 1,
//       MPI_DOUBLE,
//                 0, 8,
//                 MPI_COMM_WORLD); // правый нижний узел узел нулевого процесса

//       MPI_Recv(&Ex_local(0, Ex_local.get_Ny()),
//                static_cast<int>(Ex_local.get_Nx()), MPI_DOUBLE, 0, 9,
//                MPI_COMM_WORLD,
//                &status); // сверху (получение нулевой строки нулевого
//                процесса)
//       MPI_Recv(&Ey_local(0, Ey_local.get_Ny()),
//                static_cast<int>(Ey_local.get_Nx()), MPI_DOUBLE, 0, 10,
//                MPI_COMM_WORLD,
//                &status); // сверху (получение нулевой строки нулевого
//                процесса)
//       MPI_Recv(&Ez_local(0, Ez_local.get_Ny()),
//                static_cast<int>(Ez_local.get_Nx()), MPI_DOUBLE, 0, 11,
//                MPI_COMM_WORLD,
//                &status); // сверху (получение нулевой строки нулевого
//                процесса)

//       MPI_Recv(&Ex_local(Ex_local.get_Nx(), Ex_local.get_Ny()), 1,
//       MPI_DOUBLE,
//                0, 12, MPI_COMM_WORLD,
//                &status); // правый верхний узел (получение из нулевого
//                процесса)
//       MPI_Recv(&Ey_local(Ey_local.get_Nx(), Ey_local.get_Ny()), 1,
//       MPI_DOUBLE,
//                0, 13, MPI_COMM_WORLD,
//                &status); // правый верхний узел (получение из нулевого
//                процесса)
//       MPI_Recv(&Ez_local(Ez_local.get_Nx(), Ez_local.get_Ny()), 1,
//       MPI_DOUBLE,
//                0, 14, MPI_COMM_WORLD,
//                &status); // правый верхний узел (получение из нулевого
//                процесса)

//       MPI_Recv(&Ex_local(-1, Ex_local.get_Ny()), 1, MPI_DOUBLE, 0, 15,
//                MPI_COMM_WORLD,
//                &status); // левый верхний узел (получение из нулевого
//                процесса)
//       MPI_Recv(&Ey_local(-1, Ey_local.get_Ny()), 1, MPI_DOUBLE, 0, 16,
//                MPI_COMM_WORLD,
//                &status); // левый верхний узел (получение из нулевого
//                процесса)
//       MPI_Recv(&Ez_local(-1, Ez_local.get_Ny()), 1, MPI_DOUBLE, 0, 17,
//                MPI_COMM_WORLD,
//                &status); // левый верхний узел (получение из нулевого
//                процесса)
//     }

//     if (mpi_comm_size == 1) {
//       for (int64_t x = 0; x < Ex_local.get_Nx(); ++x) {
//         *(Ex_local.data() + x + 1) =
//             Ex_local(x, Ex_local.get_Ny() - 1);          // снизу
//         Ex_local(x, Ex_local.get_Ny()) = Ex_local(x, 0); // сверху

//         *(Ey_local.data() + x + 1) =
//             Ey_local(x, Ey_local.get_Ny() - 1);          // снизу
//         Ey_local(x, Ey_local.get_Ny()) = Ey_local(x, 0); // сверху

//         *(Ez_local.data() + x + 1) =
//             Ez_local(x, Ez_local.get_Ny() - 1);          // снизу
//         Ez_local(x, Ez_local.get_Ny()) = Ez_local(x, 0); // сверху
//       }
//       for (int64_t y = 0; y < Ex_local.get_Ny(); ++y) {
//         Ex_local(-1, y) = Ex_local(Ex_local.get_Nx() - 1, y); // слева
//         Ex_local(Ex_local.get_Nx(), y) = Ex_local(0, y);      // справа

//         Ey_local(-1, y) = Ey_local(Ey_local.get_Nx() - 1, y); // слева
//         Ey_local(Ey_local.get_Nx(), y) = Ey_local(0, y);      // справа

//         Ez_local(-1, y) = Ez_local(Ez_local.get_Nx() - 1, y); // слева
//         Ez_local(Ez_local.get_Nx(), y) = Ez_local(0, y);      // справа
//       }
//       Ex_local(-1, -1) = Ex_local(Ex_local.get_Nx() - 1,
//                                   Ex_local.get_Ny() - 1); // левый нижний
//                                   узел
//       Ex_local(-1, Ex_local.get_Ny()) =
//           Ex_local(Ex_local.get_Nx() - 1, 0); // левый верхний узел
//       Ex_local(Ex_local.get_Nx(), Ex_local.get_Ny()) =
//           Ex_local(0, 0); // правый верхний узел
//       *(Ex_local.data() + Ex_local.get_Nx() + 1) =
//           Ex_local(0, Ex_local.get_Ny() - 1); // правый нижний узел

//       Ey_local(-1, -1) = Ey_local(Ey_local.get_Nx() - 1,
//                                   Ey_local.get_Ny() - 1); // левый нижний
//                                   узел
//       Ey_local(-1, Ey_local.get_Ny()) =
//           Ey_local(Ey_local.get_Nx() - 1, 0); // левый верхний узел
//       Ey_local(Ey_local.get_Nx(), Ey_local.get_Ny()) =
//           Ey_local(0, 0); // правый верхний узел
//       *(Ey_local.data() + Ey_local.get_Nx() + 1) =
//           Ey_local(0, Ey_local.get_Ny() - 1); // правый нижний узел

//       Ez_local(-1, -1) = Ez_local(Ez_local.get_Nx() - 1,
//                                   Ez_local.get_Ny() - 1); // левый нижний
//                                   узел
//       Ez_local(-1, Ez_local.get_Ny()) =
//           Ez_local(Ez_local.get_Nx() - 1, 0); // левый верхний узел
//       Ez_local(Ez_local.get_Nx(), Ez_local.get_Ny()) =
//           Ez_local(0, 0); // правый верхний узел
//       *(Ez_local.data() + Ez_local.get_Nx() + 1) =
//           Ez_local(0, Ez_local.get_Ny() - 1); // правый нижний узел
//     }

//     // Синхронизация смежных строк (для 2-мерного случая)
//     if (mpi_comm_size > 1) {
//       if (rank != mpi_comm_size - 1 && rank != 0) {
//         MPI_Sendrecv(Ex_local.data() + Nx_local * (Ny_local - 2),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank,
//                      Ex_local.data() + Nx_local * (Ny_local - 1),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank +
//                      1, MPI_COMM_WORLD, &status); // сверху получил-отправил
//         MPI_Sendrecv(Ey_local.data() + Nx_local * (Ny_local - 2),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank,
//                      Ey_local.data() + Nx_local * (Ny_local - 1),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank +
//                      1, MPI_COMM_WORLD, &status); // сверху получил-отправил
//         MPI_Sendrecv(Ez_local.data() + Nx_local * (Ny_local - 2),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank,
//                      Ez_local.data() + Nx_local * (Ny_local - 1),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank +
//                      1, MPI_COMM_WORLD, &status); // сверху получил-отправил

//         MPI_Sendrecv(Ex_local.data() + Nx_local, static_cast<int>(Nx_local),
//                      MPI_DOUBLE, rank - 1, rank, Ex_local.data(),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank - 1, rank -
//                      1, MPI_COMM_WORLD, &status); // снизу получил-отправил
//         MPI_Sendrecv(Ey_local.data() + Nx_local, static_cast<int>(Nx_local),
//                      MPI_DOUBLE, rank - 1, rank, Ey_local.data(),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank - 1, rank -
//                      1, MPI_COMM_WORLD, &status); // снизу получил-отправил
//         MPI_Sendrecv(Ez_local.data() + Nx_local, static_cast<int>(Nx_local),
//                      MPI_DOUBLE, rank - 1, rank, Ez_local.data(),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank - 1, rank -
//                      1, MPI_COMM_WORLD, &status); // снизу получил-отправил
//       } else if (rank == 0) {
//         MPI_Sendrecv(Ex_local.data() + Nx_local * (Ny_local - 2),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank,
//                      Ex_local.data() + Nx_local * (Ny_local - 1),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank +
//                      1, MPI_COMM_WORLD, &status); // сверху получил-отправил
//         MPI_Sendrecv(Ey_local.data() + Nx_local * (Ny_local - 2),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank,
//                      Ey_local.data() + Nx_local * (Ny_local - 1),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank +
//                      1, MPI_COMM_WORLD, &status); // сверху получил-отправил
//         MPI_Sendrecv(Ez_local.data() + Nx_local * (Ny_local - 2),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank,
//                      Ez_local.data() + Nx_local * (Ny_local - 1),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank + 1, rank +
//                      1, MPI_COMM_WORLD, &status); // сверху получил-отправил
//       } else                                   // if (rank == mpi_comm_size -
//       1)
//       {
//         MPI_Sendrecv(Ex_local.data() + Nx_local, static_cast<int>(Nx_local),
//                      MPI_DOUBLE, rank - 1, rank, Ex_local.data(),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank - 1, rank -
//                      1, MPI_COMM_WORLD, &status); // снизу получил-отправил
//         MPI_Sendrecv(Ey_local.data() + Nx_local, static_cast<int>(Nx_local),
//                      MPI_DOUBLE, rank - 1, rank, Ey_local.data(),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank - 1, rank -
//                      1, MPI_COMM_WORLD, &status); // снизу получил-отправил
//         MPI_Sendrecv(Ez_local.data() + Nx_local, static_cast<int>(Nx_local),
//                      MPI_DOUBLE, rank - 1, rank, Ez_local.data(),
//                      static_cast<int>(Nx_local), MPI_DOUBLE, rank - 1, rank -
//                      1, MPI_COMM_WORLD, &status); // снизу получил-отправил
//       }
//     }

//     //
//     ==============================================================================================================================================

//     // for (i = 0; i < Nx_local; ++i)
//     //   for (j = 0; j < Ny_local; ++j)
//     //     for (k = 0; k < Nz_local; ++k)

//     for (i = 0; i < Nx_local - 2; ++i)
//       for (j = 0; j < Ny_local - 2; ++j)
//         for (k = 0; k < Nz_local; ++k) {
//           Bx_local(i, j, k) =
//               Bx_local(i, j, k) +
//               C * B_dt *
//                   (/*(Ey_local(i, j, k + 1) - Ey_local(i, j, k)) / dz*/ -(
//                        Ez_local(i, j + 1, k) - Ez_local(i, j, k)) /
//                    dy);
//           By_local(i, j, k) =
//               By_local(i, j, k) +
//               C * B_dt *
//                   ((Ez_local(i + 1, j, k) - Ez_local(i, j, k)) /
//                    dx /*- (Ex_local(i, j, k + 1) - Ex_local(i, j, k)) /
//                    dz*/);
//           Bz_local(i, j, k) =
//               Bz_local(i, j, k) +
//               C * B_dt *
//                   ((Ex_local(i, j + 1, k) - Ex_local(i, j, k)) / dy -
//                    (Ey_local(i + 1, j, k) - Ey_local(i, j, k)) / dx);
//         }
//   }

//   size_t localShiftedSize = 0;
//   if (rank == 0 || rank == mpi_comm_size - 1) {
//     localShiftedSize = Nx_local * (Ny_local - 1);
//     std::cout << "localShiftedSize(0) = " << localShiftedSize << std::endl;
//   }

//   else {
//     localShiftedSize = Nx_local * (Ny_local - 2);
//     std::cout << "localShiftedSize(drugie) = " << localShiftedSize <<
//     std::endl;
//   }
//   std::vector<double> proc_buffer_Ex(localShiftedSize);
//   std::vector<double> proc_buffer_Ey(localShiftedSize);
//   std::vector<double> proc_buffer_Ez(localShiftedSize);

//   std::vector<double> proc_buffer_Bx(localShiftedSize);
//   std::vector<double> proc_buffer_By(localShiftedSize);
//   std::vector<double> proc_buffer_Bz(localShiftedSize);

//   // std::vector<double> proc_buffer_Ex(Nx_local* (Ny_local - 1));
//   // std::vector<double> proc_buffer_Ey(Nx_local* (Ny_local - 1));
//   // std::vector<double> proc_buffer_Ez(Nx_local* (Ny_local - 1));
//   //
//   // std::vector<double> proc_buffer_Bx(Nx_local* (Ny_local - 1));
//   // std::vector<double> proc_buffer_By(Nx_local* (Ny_local - 1));
//   // std::vector<double> proc_buffer_Bz(Nx_local* (Ny_local - 1));

//   // Gather для нулевого процесса
//   if (rank == 0) {
//     // memcpy(Ex.data(), Ex_local.data(), proc_buffer_Ex.size() *
//     // sizeof(double)); memcpy(Ey.data(), Ey_local.data(),
//     proc_buffer_Ey.size()
//     // * sizeof(double)); memcpy(Ez.data(), Ez_local.data(),
//     // proc_buffer_Ez.size() * sizeof(double));

//     // memcpy(Bx.data(), Bx_local.data(), proc_buffer_Bx.size() *
//     // sizeof(double)); memcpy(By.data(), By_local.data(),
//     proc_buffer_By.size()
//     // * sizeof(double)); memcpy(Bz.data(), Bz_local.data(),
//     // proc_buffer_Bz.size() * sizeof(double));

//     memcpy(proc_buffer_Ex.data(), Ex_local.data(),
//            proc_buffer_Ex.size() * sizeof(double));
//     memcpy(proc_buffer_Ey.data(), Ey_local.data(),
//            proc_buffer_Ey.size() * sizeof(double));
//     memcpy(proc_buffer_Ez.data(), Ez_local.data(),
//            proc_buffer_Ez.size() * sizeof(double));

//     memcpy(proc_buffer_Bx.data(), Bx_local.data(),
//            proc_buffer_Bx.size() * sizeof(double));
//     memcpy(proc_buffer_By.data(), By_local.data(),
//            proc_buffer_By.size() * sizeof(double));
//     memcpy(proc_buffer_Bz.data(), Bz_local.data(),
//            proc_buffer_Bz.size() * sizeof(double));
//     // for (int i = 0; i < proc_buffer_By.size(); ++i)
//     //{
//     //   std::cout << "Bx_local[" << i << "]" << " = " << *(Bx_local.data() +
//     i)
//     //   << std::endl;
//     // }
//   }
//   // else if (rank == mpi_comm_size - 1)
//   //{
//   //   memcpy(proc_buffer_Ex.data(), Ex_local.data() + Nx_local,
//   //   proc_buffer_Ex.size() * sizeof(double)); memcpy(proc_buffer_Ey.data(),
//   //   Ey_local.data() + Nx_local, proc_buffer_Ey.size() * sizeof(double));
//   //   memcpy(proc_buffer_Ez.data(), Ez_local.data() + Nx_local,
//   //   proc_buffer_Ez.size() * sizeof(double));
//   //
//   //   memcpy(proc_buffer_Bx.data(), Bx_local.data() + Nx_local,
//   //   proc_buffer_Bx.size() * sizeof(double)); memcpy(proc_buffer_By.data(),
//   //   By_local.data() + Nx_local, proc_buffer_By.size() * sizeof(double));
//   //   memcpy(proc_buffer_Bz.data(), Bz_local.data() + Nx_local,
//   //   proc_buffer_Bz.size() * sizeof(double));
//   // }
//   else {
//     memcpy(proc_buffer_Ex.data(), Ex_local.data() + Nx_local,
//            proc_buffer_Ex.size() * sizeof(double));
//     memcpy(proc_buffer_Ey.data(), Ey_local.data() + Nx_local,
//            proc_buffer_Ey.size() * sizeof(double));
//     memcpy(proc_buffer_Ez.data(), Ez_local.data() + Nx_local,
//            proc_buffer_Ez.size() * sizeof(double));

//     memcpy(proc_buffer_Bx.data(), Bx_local.data() + Nx_local,
//            proc_buffer_Bx.size() * sizeof(double));
//     memcpy(proc_buffer_By.data(), By_local.data() + Nx_local,
//            proc_buffer_By.size() * sizeof(double));
//     memcpy(proc_buffer_Bz.data(), Bz_local.data() + Nx_local,
//            proc_buffer_Bz.size() * sizeof(double));

//     // for (int i = 0; i < proc_buffer_Bx.size(); ++i)
//     //{
//     //   std::cout << "proc_buffer_Bx[" << i << "]" << " = " <<
//     //   proc_buffer_Bx[i] << std::endl;
//     // }
//   }

//   std::vector<int> razmer{6 * 12, 5 * 12, 5 * 12};
//   std::vector<int> displs_Gather{0, 60, 96};
//   MPI_Gatherv(proc_buffer_Ex.data(), localShiftedSize, MPI_DOUBLE, Ex.data(),
//               razmer.data(), displs_Gather.data(), MPI_DOUBLE, 0,
//               MPI_COMM_WORLD);
//   MPI_Gatherv(proc_buffer_Ey.data(), localShiftedSize, MPI_DOUBLE, Ey.data(),
//               razmer.data(), displs_Gather.data(), MPI_DOUBLE, 0,
//               MPI_COMM_WORLD);
//   MPI_Gatherv(proc_buffer_Ez.data(), localShiftedSize, MPI_DOUBLE, Ez.data(),
//               razmer.data(), displs_Gather.data(), MPI_DOUBLE, 0,
//               MPI_COMM_WORLD);

//   MPI_Gatherv(proc_buffer_Bx.data(), localShiftedSize, MPI_DOUBLE, Bx.data(),
//               razmer.data(), displs_Gather.data(), MPI_DOUBLE, 0,
//               MPI_COMM_WORLD);
//   MPI_Gatherv(proc_buffer_By.data(), localShiftedSize, MPI_DOUBLE, By.data(),
//               razmer.data(), displs_Gather.data(), MPI_DOUBLE, 0,
//               MPI_COMM_WORLD);
//   MPI_Gatherv(proc_buffer_Bz.data(), localShiftedSize, MPI_DOUBLE, Bz.data(),
//               razmer.data(), displs_Gather.data(), MPI_DOUBLE, 0,
//               MPI_COMM_WORLD);

//   // MPI_Gather(proc_buffer_Ex.data(),
//   static_cast<int>(proc_buffer_Ex.size()),
//   // MPI_DOUBLE, Ex.data(), static_cast<int>(proc_buffer_Ex.size()),
//   MPI_DOUBLE,
//   // 0, MPI_COMM_WORLD); MPI_Gather(proc_buffer_Ey.data(),
//   // static_cast<int>(proc_buffer_Ey.size()), MPI_DOUBLE, Ey.data(),
//   // static_cast<int>(proc_buffer_Ey.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//   // MPI_Gather(proc_buffer_Ez.data(),
//   static_cast<int>(proc_buffer_Ez.size()),
//   // MPI_DOUBLE, Ez.data(), static_cast<int>(proc_buffer_Ez.size()),
//   MPI_DOUBLE,
//   // 0, MPI_COMM_WORLD);
//   //
//   // MPI_Gather(proc_buffer_Bx.data(),
//   static_cast<int>(proc_buffer_Bx.size()),
//   // MPI_DOUBLE, Bx.data(), static_cast<int>(proc_buffer_Bx.size()),
//   MPI_DOUBLE,
//   // 0, MPI_COMM_WORLD); MPI_Gather(proc_buffer_By.data(),
//   // static_cast<int>(proc_buffer_By.size()), MPI_DOUBLE, By.data(),
//   // static_cast<int>(proc_buffer_By.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//   // MPI_Gather(proc_buffer_Bz.data(),
//   static_cast<int>(proc_buffer_Bz.size()),
//   // MPI_DOUBLE, Bz.data(), static_cast<int>(proc_buffer_Bz.size()),
//   MPI_DOUBLE,
//   // 0, MPI_COMM_WORLD);
//   MPI_Barrier(MPI_COMM_WORLD);
//   std::cout << "We are HERE!\n" << std::endl;
//   // if (rank == 0)
//   //{
//   //   for (int i = 0; i < Bx.size(); ++i)
//   //   {
//   //     std::cout << "Bx[" << i << "]" << " = " << *(Bx.data() + i) <<
//   //     std::endl;
//   //   }
//   // }

//   if (mpi_comm_size == 1) {
//     std::memcpy(Ex.data(), Ex_local.data(),
//                 Ex_local.fullsize() * sizeof(double));
//     std::memcpy(Ey.data(), Ey_local.data(),
//                 Ey_local.fullsize() * sizeof(double));
//     std::memcpy(Ez.data(), Ez_local.data(),
//                 Ez_local.fullsize() * sizeof(double));

//     std::memcpy(Bx.data(), Bx_local.data(),
//                 Bx_local.fullsize() * sizeof(double));
//     std::memcpy(By.data(), By_local.data(),
//                 By_local.fullsize() * sizeof(double));
//     std::memcpy(Bz.data(), Bz_local.data(),
//                 Bz_local.fullsize() * sizeof(double));

//     // MPI_Send(Ex_local.data(), Nx_local * Ny_local, MPI_DOUBLE)
//   }

//   // #else
//   // for (int64_t time = 0ull; time < t; time++)
//   //{
//   //   for (i = 0; i < Nx; ++i)
//   //     for (j = 0; j < Ny; ++j)
//   //       for (k = 0; k < Nz; ++k)
//   //       {
//   //         Bx(i, j, k) = Bx(i, j, k) + C * B_dt * (/*(Ey(i, j, k + 1) -
//   Ey(i,
//   //         j, k)) / dz*/ - (Ez(i, j + 1, k) - Ez(i, j, k)) / dy); By(i, j,
//   k)
//   //         = By(i, j, k) + C * B_dt * ((Ez(i + 1, j, k) - Ez(i, j, k)) / dx
//   //         /*- (Ex(i, j, k + 1) - Ex(i, j, k)) / dz*/); Bz(i, j, k) = Bz(i,
//   j,
//   //         k) + C * B_dt * ((Ex(i, j + 1, k) - Ex(i, j, k)) / dy - (Ey(i +
//   1,
//   //         j, k) - Ey(i, j, k)) / dx);
//   //         //std::cout << "i = " << i << " j = " << j << " k = " << k <<
//   '\n';
//   //       }
//   //   // Обновляю границу - верх 2-х мерной системы
//   //   for (int64_t x = 0; x < Bx.get_Nx(); ++x)
//   //   {
//   //     *(Bx.data() + x + 1) = Bx(x, Bx.get_Ny() - 1); // снизу
//   //     Bx(x, Bx.get_Ny()) = Bx(x, 0); // сверху

//   //    *(By.data() + x + 1) = By(x, By.get_Ny() - 1); // снизу
//   //    By(x, By.get_Ny() - 1) = By(x, 0); // сверху

//   //    *(Bz.data() + x + 1) = Bz(x, Bz.get_Ny() - 1); // снизу
//   //    Bz(x, Bz.get_Ny()) = Bz(x, 0); // сверху
//   //  }
//   //  for (int64_t y = 0; y < Bx.get_Ny(); ++y)
//   //  {
//   //    Bx(-1, y) = Bx(Bx.get_Nx() - 1, y); // слева
//   //    Bx(Bx.get_Nx(), y) = Bx(0, y); // справа

//   //    By(-1, y) = By(By.get_Nx() - 1, y); // слева
//   //    By(By.get_Nx(), y) = By(0, y); // справа

//   //    Bz(-1, y) = Bz(Bz.get_Nx() - 1, y); // слева
//   //    Bz(Bz.get_Nx(), y) = Bz(0, y); // справа
//   //  }
//   //  Bx(-1, -1) = Bx(Bx.get_Nx() - 1, Bx.get_Ny() - 1); // левый нижний узел
//   //  Bx(-1, Bx.get_Ny()) = Bx(Bx.get_Nx() - 1, 0); // левый верхний узел
//   //  Bx(Bx.get_Nx(), Bx.get_Ny()) = Bx(0, 0); // правый верхний узел
//   //  *(Bx.data() + Bx.get_Nx() + 1) = Bx(0, Bx.get_Ny() - 1); // правый
//   нижний
//   //  узел

//   //  By(-1, -1) = By(By.get_Nx() - 1, By.get_Ny() - 1); // левый нижний узел
//   //  By(-1, By.get_Ny()) = By(By.get_Nx() - 1, 0); // левый верхний узел
//   //  By(By.get_Nx(), By.get_Ny()) = By(0, 0); // правый верхний узел
//   //  *(By.data() + By.get_Nx() + 1) = By(0, By.get_Ny() - 1); // правый
//   нижний
//   //  узел

//   //  Bz(-1, -1) = Bz(Bz.get_Nx() - 1, Bz.get_Ny() - 1); // левый нижний узел
//   //  Bz(-1, Bz.get_Ny()) = Bz(Bz.get_Nx() - 1, 0); // левый верхний узел
//   //  Bz(Bz.get_Nx(), Bz.get_Ny()) = Bz(0, 0); // правый верхний узел
//   //  *(Bz.data() + Bz.get_Nx() + 1) = Bz(0, Bz.get_Ny() - 1); // правый
//   нижний
//   //  узел

//   //  for (i = 0ull; i < Nx; ++i)
//   //    for (j = 0ull; j < Ny; ++j)
//   //      for (k = 0ull; k < Nz; ++k)
//   //      {
//   //        Ex(i, j, k) = Ex(i, j, k) + C * E_dt * ((Bz(i, j, k) - Bz(i, j -
//   1,
//   //        k)) / dy /*- (By(i, j, k) - By(i, j, k - 1)) / dz*/); Ey(i, j, k)
//   =
//   //        Ey(i, j, k) + C * E_dt * (/*(Bx(i, j, k) - Bx(i, j, k - 1)) /
//   dz*/ -
//   //        (Bz(i, j, k) - Bz(i - 1, j, k)) / dx); Ez(i, j, k) = Ez(i, j, k)
//   + C
//   //        * E_dt * ((By(i, j, k) - By(i - 1, j, k)) / dx - (Bx(i, j, k) -
//   //        Bx(i, j - 1, k)) / dy);
//   //      }
//   //  // Обновляю границу - верх 2-х мерной системы
//   //  for (int64_t x = 0; x < Ex.get_Nx(); ++x)
//   //  {
//   //    *(Ex.data() + x + 1) = Ex(x, Ex.get_Ny() - 1); // снизу
//   //    Ex(x, Ex.get_Ny()) = Ex(x, 0); // сверху

//   //    *(Ey.data() + x + 1) = Ey(x, Ey.get_Ny() - 1); // снизу
//   //    Ey(x, Ey.get_Ny()) = Ey(x, 0); // сверху

//   //    *(Ez.data() + x + 1) = Ez(x, Ez.get_Ny() - 1); // снизу
//   //    Ez(x, Ez.get_Ny()) = Ez(x, 0); // сверху
//   //  }
//   //  for (int64_t y = 0; y < Ex.get_Ny(); ++y)
//   //  {
//   //    Ex(-1, y) = Ex(Ex.get_Nx() - 1, y); // слева
//   //    Ex(Ex.get_Nx(), y) = Ex(0, y); // справа

//   //    Ey(-1, y) = Ey(Ey.get_Nx() - 1, y); // слева
//   //    Ey(Ey.get_Nx(), y) = Ey(0, y); // справа

//   //    Ez(-1, y) = Ez(Ez.get_Nx() - 1, y); // слева
//   //    Ez(Ez.get_Nx(), y) = Ez(0, y); // справа
//   //  }
//   //  Ex(-1, -1) = Ex(Ex.get_Nx() - 1, Ex.get_Ny() - 1); // левый нижний узел
//   //  Ex(-1, Ex.get_Ny()) = Ex(Ex.get_Nx() - 1, 0); // левый верхний узел
//   //  Ex(Ex.get_Nx(), Ex.get_Ny()) = Ex(0, 0); // правый верхний узел
//   //  *(Ex.data() + Ex.get_Nx() + 1) = Ex(0, Ex.get_Ny() - 1); // правый
//   нижний
//   //  узел

//   //  Ey(-1, -1) = Ey(Ey.get_Nx() - 1, Ey.get_Ny() - 1); // левый нижний узел
//   //  Ey(-1, Ey.get_Ny()) = Ey(Ey.get_Nx() - 1, 0); // левый верхний узел
//   //  Ey(Ey.get_Nx(), Ey.get_Ny()) = Ey(0, 0); // правый верхний узел
//   //  *(Ey.data() + Ey.get_Nx() + 1) = Ey(0, Ey.get_Ny() - 1); // правый
//   нижний
//   //  узел

//   //  Ez(-1, -1) = Ez(Ez.get_Nx() - 1, Ez.get_Ny() - 1); // левый нижний узел
//   //  Ez(-1, Ez.get_Ny()) = Ez(Ez.get_Nx() - 1, 0); // левый верхний узел
//   //  Ez(Ez.get_Nx(), Ez.get_Ny()) = Ez(0, 0); // правый верхний узел
//   //  *(Ez.data() + Ez.get_Nx() + 1) = Ez(0, Ez.get_Ny() - 1); // правый
//   нижний
//   //  узел

//   //  for (i = 0ull; i < Nx; ++i)
//   //    for (j = 0ull; j < Ny; ++j)
//   //      for (k = 0ull; k < Nz; ++k)
//   //      {
//   //        Bx(i, j, k) = Bx(i, j, k) + C * B_dt * (/*(Ey(i, j, k + 1) -
//   Ey(i,
//   //        j, k)) / dz*/ - (Ez(i, j + 1, k) - Ez(i, j, k)) / dy); By(i, j,
//   k) =
//   //        By(i, j, k) + C * B_dt * ((Ez(i + 1, j, k) - Ez(i, j, k)) / dx
//   /*-
//   //        (Ex(i, j, k + 1) - Ex(i, j, k)) / dz*/); Bz(i, j, k) = Bz(i, j,
//   k) +
//   //        C * B_dt * ((Ex(i, j + 1, k) - Ex(i, j, k)) / dy - (Ey(i + 1, j,
//   k)
//   //        - Ey(i, j, k)) / dx);
//   //      }
//   //  // Обновляю границу - верх 2-х мерной системы
//   //  for (int64_t x = 0; x < Bx.get_Nx(); ++x)
//   //  {
//   //    *(Bx.data() + x + 1) = Bx(x, Bx.get_Ny() - 1); // снизу
//   //    Bx(x, Bx.get_Ny()) = Bx(x, 0); // сверху

//   //    *(By.data() + x + 1) = By(x, By.get_Ny() - 1); // снизу
//   //    By(x, By.get_Ny()) = By(x, 0); // сверху

//   //    *(Bz.data() + x + 1) = Bz(x, Bz.get_Ny() - 1); // снизу
//   //    Bz(x, Bz.get_Ny()) = Bz(x, 0); // сверху
//   //  }
//   //  for (int64_t y = 0; y < Bx.get_Ny(); ++y)
//   //  {
//   //    Bx(-1, y) = Bx(Bx.get_Nx() - 1, y); // слева
//   //    Bx(Bx.get_Nx(), y) = Bx(0, y); // справа

//   //    By(-1, y) = By(By.get_Nx() - 1, y); // слева
//   //    By(By.get_Nx(), y) = By(0, y); // справа

//   //    Bz(-1, y) = Bz(Bz.get_Nx() - 1, y); // слева
//   //    Bz(Bz.get_Nx(), y) = Bz(0, y); // справа
//   //  }
//   //  Bx(-1, -1) = Bx(Bx.get_Nx() - 1, Bx.get_Ny() - 1); // левый нижний узел
//   //  Bx(-1, Bx.get_Ny()) = Bx(Bx.get_Nx() - 1, 0); // левый верхний узел
//   //  Bx(Bx.get_Nx(), Bx.get_Ny()) = Bx(0, 0); // правый верхний узел
//   //  *(Bx.data() + Bx.get_Nx() + 1) = Bx(0, Bx.get_Ny() - 1); // правый
//   нижний
//   //  узел

//   //  By(-1, -1) = By(By.get_Nx() - 1, By.get_Ny() - 1); // левый нижний узел
//   //  By(-1, By.get_Ny()) = By(By.get_Nx() - 1, 0); // левый верхний узел
//   //  By(By.get_Nx(), By.get_Ny()) = By(0, 0); // правый верхний узел
//   //  *(By.data() + By.get_Nx() + 1) = By(0, By.get_Ny() - 1); // правый
//   нижний
//   //  узел

//   //  Bz(-1, -1) = Bz(Bz.get_Nx() - 1, Bz.get_Ny() - 1); // левый нижний узел
//   //  Bz(-1, Bz.get_Ny()) = Bz(Bz.get_Nx() - 1, 0); // левый верхний узел
//   //  Bz(Bz.get_Nx(), Bz.get_Ny()) = Bz(0, 0); // правый верхний узел
//   //  *(Bz.data() + Bz.get_Nx() + 1) = Bz(0, Bz.get_Ny() - 1); // правый
//   нижний
//   //  узел
//   //}

//   // #endif // MPI
// }

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
    std::cout << "\nGet_Axis: Error! Wrong Components!\n";
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
    std::cout << "\nGet_Axis: Error! Wrong Components!\n";
    exit(-1);
  }
}

void FDTD::FDTD::boundary_synchronization(/*Сюда добавить поля*/) {

  // СИСТЕМА ТЕГОВ ДЛЯ MPI (ОТПРАВКА):
  // 0 - низ
  // 1 - верх

  int rank = 0;
  int mpi_comm_size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_comm_size);

  if (rank == 0) {
    // std::cout << "Nx = " << Nx << " Ny = " << Ny << " Nz = " << Nz << std::endl;
  }


  // Сначала рассматривается разбиение по OY для 2D

  // Слева и справа (граничные условия для 2D) (Всё это внутри процесса, поэтому без MPI сообщений)
  for (int64_t i = 0; i < Ny; ++i) {
    // Слева
    Ex(-1, i) = Ex(Nx - 1, i);
    Ey(-1, i) = Ey(Nx - 1, i);
    Ez(-1, i) = Ez(Nx - 1, i);
    Bx(-1, i) = Bx(Nx - 1, i);
    By(-1, i) = By(Nx - 1, i);
    Bz(-1, i) = Bz(Nx - 1, i);

    // Справа
    Ex(Nx, i) = Ex(0, i);
    Ey(Nx, i) = Ey(0, i);
    Ez(Nx, i) = Ez(0, i);
    Bx(Nx, i) = Bx(0, i);
    By(Nx, i) = By(0, i);
    Bz(Nx, i) = Bz(0, i);
  }

  // В NY фактически хранится Ny_local (то есть размер без оболочки)

  // Нужно потом будет учесть ситуацию, когда только один процесс
  // ТАКЖЕ ИНДЕКСЫ ЗДЕСЬ НЕ УЧИТЫВАЮТ OZ НАПРАВЛЕНИЕ!
  // НЕ забыть добавить правильную обработку угловых значений для 0 и последнего процессов
  // ****************************************************************************************************
  // Низ каждого процесса передаётся вверх предыдущему процессу (P.S. не нулевая строка, а первая)

  if (mpi_comm_size > 1) {

    if (rank != 0) {
      MPI_Send(&Ex(-1, 0), Nx + 2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
      MPI_Send(&Ey(-1, 0), Nx + 2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
      MPI_Send(&Ez(-1, 0), Nx + 2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
      MPI_Send(&Bx(-1, 0), Nx + 2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
      MPI_Send(&By(-1, 0), Nx + 2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
      MPI_Send(&Bz(-1, 0), Nx + 2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
    }

    // Принятие верха от низа следующего процесса (P.S. последняя строка)
    if (rank != mpi_comm_size - 1) {
      MPI_Recv(&Ex(-1, Ny), Nx + 2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&Ey(-1, Ny), Nx + 2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&Ez(-1, Ny), Nx + 2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&Bx(-1, Ny), Nx + 2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&By(-1, Ny), Nx + 2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&Bz(-1, Ny), Nx + 2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // ****************************************************************************************************

    // Верх каждого процесса передаётся вниз следующему процессу (P.S. не последняя строка, а предпоследняя)
    if (rank != mpi_comm_size - 1) {
      MPI_Send(&Ex(-1, Ny - 1), Nx + 2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
      MPI_Send(&Ey(-1, Ny - 1), Nx + 2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
      MPI_Send(&Ez(-1, Ny - 1), Nx + 2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
      MPI_Send(&Bx(-1, Ny - 1), Nx + 2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
      MPI_Send(&By(-1, Ny - 1), Nx + 2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
      MPI_Send(&Bz(-1, Ny - 1), Nx + 2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
    }
    
    // Принятие низа от верха предыдущего процесса (P.S нулевая строка)
    if (rank != 0) {
      MPI_Recv(Ex.data(), Nx + 2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(Ey.data(), Nx + 2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(Ez.data(), Nx + 2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(Bx.data(), Nx + 2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(By.data(), Nx + 2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(Bz.data(), Nx + 2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // ****************************************************************************************************

    // Синхронизация для нулевого и последнего процессов

    // Отправка низа вверх последнего процесса (P.S. не нулевая строка, а первая), а также получение низа (нулевой строки) от последнего процесса
    if (rank == 0) {
      MPI_Sendrecv(&Ex(-1, 0), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 0, Ex.data(), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Sendrecv(&Ey(-1, 0), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 0, Ey.data(), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Sendrecv(&Ez(-1, 0), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 0, Ez.data(), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Sendrecv(&Bx(-1, 0), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 0, Bx.data(), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Sendrecv(&By(-1, 0), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 0, By.data(), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Sendrecv(&Bz(-1, 0), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 0, Bz.data(), Nx + 2, MPI_DOUBLE, mpi_comm_size - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Отправка верха вниз нулевого процесса (P.S. не последняя строка, а предпоследняя), а также получение верха (последней строки) от нулевого процесса
    if (rank == mpi_comm_size - 1) {
      MPI_Sendrecv(&Ex(-1, Ny - 1), Nx + 2, MPI_DOUBLE, 0, 1, &Ex(-1, Ny), Nx + 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Sendrecv(&Ey(-1, Ny - 1), Nx + 2, MPI_DOUBLE, 0, 1, &Ey(-1, Ny), Nx + 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Sendrecv(&Ez(-1, Ny - 1), Nx + 2, MPI_DOUBLE, 0, 1, &Ez(-1, Ny), Nx + 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Sendrecv(&Bx(-1, Ny - 1), Nx + 2, MPI_DOUBLE, 0, 1, &Bx(-1, Ny), Nx + 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Sendrecv(&By(-1, Ny - 1), Nx + 2, MPI_DOUBLE, 0, 1, &By(-1, Ny), Nx + 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Sendrecv(&Bz(-1, Ny - 1), Nx + 2, MPI_DOUBLE, 0, 1, &Bz(-1, Ny), Nx + 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  } else {
    // Верх (с угловыми узлами)
    for (int64_t x = -1; x < Nx + 1; ++x) {
      Ex(x, Ny) = Ex(x, 0);
      Ey(x, Ny) = Ey(x, 0);
      Ez(x, Ny) = Ez(x, 0);
      Bx(x, Ny) = Bx(x, 0);
      By(x, Ny) = By(x, 0);
      Bz(x, Ny) = Bz(x, 0);
    }
    // Низ (с угловыми узлами)
    for (int64_t x = -1; x < Nx + 1; ++x) {
      Ex(x, -1) = Ex(x, Ny - 1);
      Ey(x, -1) = Ey(x, Ny - 1);
      Ez(x, -1) = Ez(x, Ny - 1);
      Bx(x, -1) = Bx(x, Ny - 1);
      By(x, -1) = By(x, Ny - 1);
      Bz(x, -1) = Bz(x, Ny - 1);
    }
  }


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

void FDTD::FDTD::boundary_synchronization_3D() {

  // Для 3д случая добавлен цикл по z для синхронизации каждого яруса 2д границ \
  а также отдельно без цикла добавлена синхронизация нижней и верхней граней \
  относительно Oz.

  // Для самого нижнего -1 и верхнего Nz + 1 яруса, необходимо сначала сам центр
  // синхронизировать, а затем уже боковые узлы, так как они зависят от центра.

  // ========================================================================

  // ОБНОВЛЕНИЕ НИЖНЕЙ ГРАНИЦЫ ОТНОСИТЕЛЬНО Oz
   // Отдельно угловые узлы в этом случае не обрабатываю, так как \
если обновить бока (слева и справа), то для угловых узлов будет всё включено

  // Центр нижней грани
#pragma omp for collapse(2)
  for (int64_t x = 0; x < Nx; ++x) {
    for (int64_t y = 0; y < Ny; ++y) {
      Ex(x, y, -1) = Ex(x, y, Nz - 1);
      Ey(x, y, -1) = Ey(x, y, Nz - 1);
      Ez(x, y, -1) = Ez(x, y, Nz - 1);
      Bx(x, y, -1) = Bx(x, y, Nz - 1);
      By(x, y, -1) = By(x, y, Nz - 1);
      Bz(x, y, -1) = Bz(x, y, Nz - 1);
    }
  }
  // Здесь вот не на 100% уверен насчёт Nz - 1, но вроде так
  // Слева нижней грани (за исключением угловых узлов)
#pragma omp for
  for (int64_t y = 0; y < Ny; ++y) {
    Ex(-1, y, -1) = Ex(Nx - 1, y, Nz - 1);
    Ey(-1, y, -1) = Ey(Nx - 1, y, Nz - 1);
    Ez(-1, y, -1) = Ez(Nx - 1, y, Nz - 1);
    Bx(-1, y, -1) = Bx(Nx - 1, y, Nz - 1);
    By(-1, y, -1) = By(Nx - 1, y, Nz - 1);
    Bz(-1, y, -1) = Bz(Nx - 1, y, Nz - 1);
  }
  // Справа нижней грани (за исключением угловых узлов)
#pragma omp for
  for (int64_t y = 0; y < Ny; ++y) {
    Ex(Nx, y, -1) = Ex(0, y, Nz - 1);
    Ey(Nx, y, -1) = Ey(0, y, Nz - 1);
    Ez(Nx, y, -1) = Ez(0, y, Nz - 1);
    Bx(Nx, y, -1) = Bx(0, y, Nz - 1);
    By(Nx, y, -1) = By(0, y, Nz - 1);
    Bz(Nx, y, -1) = Bz(0, y, Nz - 1);
  }
  // Верх нижней грани (с учётом угловых узлов)
#pragma omp for
  for (int64_t x = -1; x < Nx + 1; ++x) {
    Ex(x, Ny, -1) = Ex(x, 0, Nz - 1);
    Ey(x, Ny, -1) = Ey(x, 0, Nz - 1);
    Ez(x, Ny, -1) = Ez(x, 0, Nz - 1);
    Bx(x, Ny, -1) = Bx(x, 0, Nz - 1);
    By(x, Ny, -1) = By(x, 0, Nz - 1);
    Bz(x, Ny, -1) = Bz(x, 0, Nz - 1);
  }
  // Низ нижней грани (с учётом угловых узлов)
#pragma omp for
  for (int64_t x = -1; x < Nx + 1; ++x) {
    Ex(x, -1, -1) = Ex(x, Ny - 1, Nz - 1);
    Ey(x, -1, -1) = Ey(x, Ny - 1, Nz - 1);
    Ez(x, -1, -1) = Ez(x, Ny - 1, Nz - 1);
    Bx(x, -1, -1) = Bx(x, Ny - 1, Nz - 1);
    By(x, -1, -1) = By(x, Ny - 1, Nz - 1);
    Bz(x, -1, -1) = Bz(x, Ny - 1, Nz - 1);
  }
  // ========================================================================
  // ОБНОВЛЕНИЕ ВЕРХНЕЙ ГРАНИЦЫ ОТНОСИТЕЛЬНО Oz

  // Центр верхней грани
#pragma omp for collapse(2)
  for (int64_t x = 0; x < Nx; ++x) {
    for (int64_t y = 0; y < Ny; ++y) {
      Ex(x, y, Nz) = Ex(x, y, 0);
      Ey(x, y, Nz) = Ey(x, y, 0);
      Ez(x, y, Nz) = Ez(x, y, 0);
      Bx(x, y, Nz) = Bx(x, y, 0);
      By(x, y, Nz) = By(x, y, 0);
      Bz(x, y, Nz) = Bz(x, y, 0);
    }
  }
  // Слева верхней грани (за исключением угловых узлов)
#pragma omp for
  for (int64_t y = 0; y < Ny; ++y) {
    Ex(-1, y, Nz) = Ex(Nx - 1, y, 0);
    Ey(-1, y, Nz) = Ey(Nx - 1, y, 0);
    Ez(-1, y, Nz) = Ez(Nx - 1, y, 0);
    Bx(-1, y, Nz) = Bx(Nx - 1, y, 0);
    By(-1, y, Nz) = By(Nx - 1, y, 0);
    Bz(-1, y, Nz) = Bz(Nx - 1, y, 0);
  }
  // Справа верхней грани (за исключением угловых узлов)
#pragma omp for
  for (int64_t y = 0; y < Ny; ++y) {
    Ex(Nx, y, Nz) = Ex(0, y, 0);
    Ey(Nx, y, Nz) = Ey(0, y, 0);
    Ez(Nx, y, Nz) = Ez(0, y, 0);
    Bx(Nx, y, Nz) = Bx(0, y, 0);
    By(Nx, y, Nz) = By(0, y, 0);
    Bz(Nx, y, Nz) = Bz(0, y, 0);
  }
  // Верх верхней грани (с учётом угловых узлов)
#pragma omp for
  for (int64_t x = -1; x < Nx + 1; ++x) {
    Ex(x, Ny, Nz) = Ex(x, 0, 0);
    Ey(x, Ny, Nz) = Ey(x, 0, 0);
    Ez(x, Ny, Nz) = Ez(x, 0, 0);
    Bx(x, Ny, Nz) = Bx(x, 0, 0);
    By(x, Ny, Nz) = By(x, 0, 0);
    Bz(x, Ny, Nz) = Bz(x, 0, 0);
  }
  // Низ верхней грани (с учётом угловых узлов)
#pragma omp for
  for (int64_t x = -1; x < Nx + 1; ++x) {
    Ex(x, -1, Nz) = Ex(x, Ny - 1, 0);
    Ey(x, -1, Nz) = Ey(x, Ny - 1, 0);
    Ez(x, -1, Nz) = Ez(x, Ny - 1, 0);
    Bx(x, -1, Nz) = Bx(x, Ny - 1, 0);
    By(x, -1, Nz) = By(x, Ny - 1, 0);
    Bz(x, -1, Nz) = Bz(x, Ny - 1, 0);
  }

  // ========================================================================
  // Синхронизация границ для 2д случая
  // for (int64_t z{-1}; z < Nz + 1; ++z) { странно. Вроде как должно быть 0, Nz
  // Щас попробую так

  // for (int64_t z{-1}; z < Nz + 1; ++z) {
#pragma omp for
  for (int64_t z = 0; z < Nz; ++z) {
    // Слева 2D (за исключением угловых узлов)
    for (int64_t y = 0; y < Ny; ++y) {
      Ex(-1, y, z) = Ex(Nx - 1, y, z);
      Ey(-1, y, z) = Ey(Nx - 1, y, z);
      Ez(-1, y, z) = Ez(Nx - 1, y, z);
      Bx(-1, y, z) = Bx(Nx - 1, y, z);
      By(-1, y, z) = By(Nx - 1, y, z);
      Bz(-1, y, z) = Bz(Nx - 1, y, z);
    }
    // Справа 2D (за исключением угловых узлов)
    for (int64_t y = 0; y < Ny; ++y) {
      Ex(Nx, y, z) = Ex(0, y, z);
      Ey(Nx, y, z) = Ey(0, y, z);
      Ez(Nx, y, z) = Ez(0, y, z);
      Bx(Nx, y, z) = Bx(0, y, z);
      By(Nx, y, z) = By(0, y, z);
      Bz(Nx, y, z) = Bz(0, y, z);
    }
    // Верх 2D (с учётом угловых узлов)
    for (int64_t x = -1; x < Nx + 1; ++x) {
      Ex(x, Ny, z) = Ex(x, 0, z);
      Ey(x, Ny, z) = Ey(x, 0, z);
      Ez(x, Ny, z) = Ez(x, 0, z);
      Bx(x, Ny, z) = Bx(x, 0, z);
      By(x, Ny, z) = By(x, 0, z);
      Bz(x, Ny, z) = Bz(x, 0, z);
    }
    // Низ 2D (с учётом угловых узлов)
    for (int64_t x = -1; x < Nx + 1; ++x) {
      Ex(x, -1, z) = Ex(x, Ny - 1, z);
      Ey(x, -1, z) = Ey(x, Ny - 1, z);
      Ez(x, -1, z) = Ez(x, Ny - 1, z);
      Bx(x, -1, z) = Bx(x, Ny - 1, z);
      By(x, -1, z) = By(x, Ny - 1, z);
      Bz(x, -1, z) = Bz(x, Ny - 1, z);
    }
  }
}
