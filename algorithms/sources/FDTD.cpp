#include "FDTD.hpp"

FDTD::FDTD::FDTD(const std::tuple<int64_t, int64_t, int64_t>& Nx_Ny_Nz,
  const std::tuple<double, double, double>& ax_ay_az,
  const std::tuple<double, double, double>& bx_by_bz, double _dt)
  : Ex(std::get<0>(Nx_Ny_Nz), std::get<1>(Nx_Ny_Nz), std::get<2>(Nx_Ny_Nz)),
  Ey(std::get<0>(Nx_Ny_Nz), std::get<1>(Nx_Ny_Nz), std::get<2>(Nx_Ny_Nz)),
  Ez(std::get<0>(Nx_Ny_Nz), std::get<1>(Nx_Ny_Nz), std::get<2>(Nx_Ny_Nz)),
  Bx(std::get<0>(Nx_Ny_Nz), std::get<1>(Nx_Ny_Nz), std::get<2>(Nx_Ny_Nz)),
  By(std::get<0>(Nx_Ny_Nz), std::get<1>(Nx_Ny_Nz), std::get<2>(Nx_Ny_Nz)),
  Bz(std::get<0>(Nx_Ny_Nz), std::get<1>(Nx_Ny_Nz), std::get<2>(Nx_Ny_Nz)) {
  
  Nx = std::get<0>(Nx_Ny_Nz);
  Ny = std::get<1>(Nx_Ny_Nz);
  Nz = std::get<2>(Nx_Ny_Nz);

  ax = std::get<0>(ax_ay_az);
  ay = std::get<1>(ax_ay_az);
  az = std::get<2>(ax_ay_az);

  bx = std::get<0>(bx_by_bz);
  by = std::get<1>(bx_by_bz);
  bz = std::get<2>(bx_by_bz);

  dx = (bx - ax) / Nx;
  dy = (by - ay) / Ny;
  dz = (bz - az) / Nz;

  dt = _dt;
}


FDTD::FDTD::FDTD(const FDTD& _fields) 
  : Ex(_fields.Ex),
  Ey(_fields.Ey),
  Ez(_fields.Ez),
  Bx(_fields.Bx),
  By(_fields.By),
  Bz(_fields.Bz) {

  Nx = _fields.Nx;
  Ny = _fields.Ny;
  Nz = _fields.Nz;

  ax = _fields.ax;
  ay = _fields.ay;
  az = _fields.az;

  bx = _fields.bx;
  by = _fields.by;
  bz = _fields.bz;

  dx = _fields.dx;
  dy = _fields.dy;
  dz = _fields.dz;

  dt = _fields.dt;
}

FDTD::FDTD::FDTD(FDTD&& _fields) noexcept
  : Ex(std::move(_fields.Ex)),
  Ey(std::move(_fields.Ey)),
  Ez(std::move(_fields.Ez)),
  Bx(std::move(_fields.Bx)),
  By(std::move(_fields.By)),
  Bz(std::move(_fields.Bz)) {

  Nx = _fields.Nx;
  Ny = _fields.Ny;
  Nz = _fields.Nz;

  ax = _fields.ax;
  ay = _fields.ay;
  az = _fields.az;

  bx = _fields.bx;
  by = _fields.by;
  bz = _fields.bz;

  dx = _fields.dx;
  dy = _fields.dy;
  dz = _fields.dz;

  dt = _fields.dt;
}

FDTD::FDTD& FDTD::FDTD::operator=(const FDTD& _fields)
{
  if (this != &_fields)
  {
    this->FDTD::FDTD(_fields);
  }
  return *this;
}

FDTD::FDTD& FDTD::FDTD::operator=(FDTD&& _fields) noexcept
{
  if (this != &_fields)
  {
    this->FDTD::FDTD(std::move(_fields));
  }
  return *this;
}

void FDTD::FDTD::set_Nx_Ny_Nz(int64_t _Nx, int64_t _Ny, int64_t _Nz)
{
  Nx = _Nx;
  Ny = _Ny;

  if (!Nz) _Nz = 1ull;
  else Nz = _Nz;

  dx = (bx - ax) / Nx;
  dy = (by - ay) / Ny;
  dz = (bz - az) / Nz;

  Ex.resize_field(_Nx, _Ny, _Nz);
  Ey.resize_field(_Nx, _Ny, _Nz);
  Ez.resize_field(_Nx, _Ny, _Nz);
  Bx.resize_field(_Nx, _Ny, _Nz);
  By.resize_field(_Nx, _Ny, _Nz);
  Bz.resize_field(_Nx, _Ny, _Nz);
}

void FDTD::FDTD::field_update(const double t) {
  if (dt == 0.0)
  {
    std::cout << "Time step is null";
    exit(-1);
  }

  int64_t i = 0ull;
  int64_t j = 0ull;

  for (double time = 0.0; time < t; time += dt)
  {
    for (j = 0ull; j < Ny; ++j)
      for (i = 0ull; i < Nx; ++i) {
        Ex(i, j) =
          Ex(i, j) + C * dt * 0.5 * ((Bz(i, j + 1ull) - Bz(i, j - 1ull)) / dy);
        Ey(i, j) =
          Ey(i, j) - C * dt * 0.5 * ((Bz(i + 1ull, j) - Bz(i - 1ull, j)) / dx);
        Ez(i, j) = Ez(i, j) + C * dt * 0.5 *
          (((By(i + 1ull, j) - By(i - 1ull, j)) / dx) - (Bx(i, j + 1ull) - Bx(i, j - 1ull)) / dy);
      }
    for (j = 0ull; j < Ny; ++j)
      for (i = 0ull; i < Nx; ++i) {
        Bx(i, j) =
          Bx(i, j) - C * dt * 0.5 * ((Ez(i, j + 1ull) - Ez(i, j - 1ull)) / dy);
        By(i, j) =
          By(i, j) + C * dt * 0.5 * ((Ez(i + 1ull, j) - Ez(i - 1ull, j)) / dx);
        Bz(i, j) = Bz(i, j) - C * dt * 0.5 *
          (((Ey(i + 1ull, j) - Ey(i - 1ull, j)) / dx) -
            (Ex(i, j + 1ull) - Ex(i, j - 1ull)) / dy);
      }
  }
  //for (double time = 0.0; time < t; time += dt)
  //  for (int64_t j = 0; j < Ny; ++j)
  //    for (int64_t i = 0; i < Nx; ++i)
  //    {
  //      Ex(i, j) = Ex(i, j) + C * dt * ((Bz(i, j + 1) - Bz(i, j - 1)) / (2.0 * dy));
  //      Ey(i, j) = Ey(i, j) - C * dt * ((Bz(i + 1, j) - Bz(i - 1, j)) / (2.0 * dy));
  //      Ez(i, j) = Ez(i, j) + C * dt * (((By(i + 1, j) - By(i - 1, j)) / (2.0 * dx)) - ((Bx(i, j + 1) - (Bx(i, j - 1))) / (2.0 * dy)));

  //      Bx(i, j) = Bx(i, j) - C * dt * ((Ez(i, j + 1) - (Ez(i, j - 1))) / (2.0 * dy));
  //      By(i, j) = By(i, j) + C * dt * ((Ez(i + 1, j) - (Ez(i - 1, j))) / (2.0 * dx));
  //      Bz(i, j) = Bz(i, j) - C * dt * (((Ey(i + 1, j) - Ey(i - 1, j)) / (2.0 * dx)) - ((Ex(i, j + 1) - (Ex(i, j - 1))) / (2.0 * dy)));
  //    }
}

void FDTD::FDTD::field_update(const int64_t t)
{
  if (dt == 0.0)
  {
    std::cout << "Time step is null";
    exit(-1);
  }

  int64_t i = 0ull;
  int64_t j = 0ull;

  for (int64_t time = 0ull; time < t; time++)
  {
    for (j = 0ull; j < Ny; ++j)
      for (i = 0ull; i < Nx; ++i) {
        Ex(i, j) =
          Ex(i, j) + C * dt * 0.5 * ((Bz(i, j + 1ull) - Bz(i, j - 1ull)) / dy);
        Ey(i, j) =
          Ey(i, j) - C * dt * 0.5 * ((Bz(i + 1ull, j) - Bz(i - 1ull, j)) / dx);
        Ez(i, j) = Ez(i, j) + C * dt * 0.5 *
          (((By(i + 1ull, j) - By(i - 1ull, j)) / dx) - (Bx(i, j + 1ull) - Bx(i, j - 1ull)) / dy);
      }
    for (j = 0ull; j < Ny; ++j)
      for (i = 0ull; i < Nx; ++i) {
        Bx(i, j) =
          Bx(i, j) - C * dt * 0.5 * ((Ez(i, j + 1ull) - Ez(i, j - 1ull)) / dy);
        By(i, j) =
          By(i, j) + C * dt * 0.5 * ((Ez(i + 1ull, j) - Ez(i - 1ull, j)) / dx);
        Bz(i, j) = Bz(i, j) - C * dt * 0.5 *
          (((Ey(i + 1ull, j) - Ey(i - 1ull, j)) / dx) -
            (Ex(i, j + 1ull) - Ex(i, j - 1ull)) / dy);
      }
  }
}

void FDTD::FDTD::shifted_field_update(const double t)
{
  if (dt == 0.0)
  {
    std::cout << "Time step is null";
    exit(-1);
  }

  double B_dt = dt * 0.5;
  double E_dt = dt;
  int64_t i = 0ull;
  int64_t j = 0ull;
  int64_t k = 0ull;

  //double start = omp_get_wtime();
  for (double time = 0.0; time < t; time += dt) // ПОПРАВИТЬ НА time += E_dt
  {
    for (i = 0ull; i < Nx; ++i)
      for (j = 0ull; j < Ny; ++j)
        for (k = 0ull; k < Nz; ++k)
        {
          Bx(i, j, k) = Bx(i, j, k) + C * B_dt * ((Ey(i, j, k + 1) - Ey(i, j, k)) / dz - (Ez(i, j + 1, k) - Ez(i, j, k)) / dy);
          By(i, j, k) = By(i, j, k) + C * B_dt * ((Ez(i + 1, j, k) - Ez(i, j, k)) / dx - (Ex(i, j, k + 1) - Ex(i, j, k)) / dz);
          Bz(i, j, k) = Bz(i, j, k) + C * B_dt * ((Ex(i, j + 1, k) - Ex(i, j, k)) / dy - (Ey(i + 1, j, k) - Ey(i, j, k)) / dx);
        }
    for (i = 0ull; i < Nx; ++i)
      for (j = 0ull; j < Ny; ++j)
        for (k = 0ull; k < Nz; ++k)
        {
          Ex(i, j, k) = Ex(i, j, k) + C * E_dt * ((Bz(i, j, k) - Bz(i, j - 1, k)) / dy - (By(i, j, k) - By(i, j, k - 1)) / dz);
          Ey(i, j, k) = Ey(i, j, k) + C * E_dt * ((Bx(i, j, k) - Bx(i, j, k - 1)) / dz - (Bz(i, j, k) - Bz(i - 1, j, k)) / dx);
          Ez(i, j, k) = Ez(i, j, k) + C * E_dt * ((By(i, j, k) - By(i - 1, j, k)) / dx - (Bx(i, j, k) - Bx(i, j - 1, k)) / dy);
        }
    for (i = 0ull; i < Nx; ++i)
      for (j = 0ull; j < Ny; ++j)
        for (k = 0ull; k < Nz; ++k)
        {
          Bx(i, j, k) = Bx(i, j, k) + C * B_dt * ((Ey(i, j, k + 1) - Ey(i, j, k)) / dz - (Ez(i, j + 1, k) - Ez(i, j, k)) / dy);
          By(i, j, k) = By(i, j, k) + C * B_dt * ((Ez(i + 1, j, k) - Ez(i, j, k)) / dx - (Ex(i, j, k + 1) - Ex(i, j, k)) / dz);
          Bz(i, j, k) = Bz(i, j, k) + C * B_dt * ((Ex(i, j + 1, k) - Ex(i, j, k)) / dy - (Ey(i + 1, j, k) - Ey(i, j, k)) / dx);
        }
  }
  //double end = omp_get_wtime();
  //std::cout << "Time = " << end - start << std::endl;
}

void FDTD::FDTD::shifted_field_update(const int64_t t)
{
  if (dt == 0.0)
  {
    std::cout << "Time step is null";
    exit(-1);
  }

  double B_dt = dt * 0.5;
  double E_dt = dt;
  int64_t i = 0ull;
  int64_t j = 0ull;
  int64_t k = 0ull;
  //double start = omp_get_wtime();

//#ifdef MPI

  int rank = 0, mpi_comm_size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_comm_size);
  MPI_Status status;

  std::cout << "MPI size = " << mpi_comm_size << '\n';

  int64_t Nx_local = get_Nx() + 2;
  std::cout << "Eto NY = " << get_Ny() << std::endl;
  int64_t Ny_local{};
  if (mpi_comm_size == 1) 
    Ny_local = (get_Ny() + 2) / mpi_comm_size;
  else 
    Ny_local = ((get_Ny() + 2) / mpi_comm_size) + 1; //   int64_t Ny_local = (get_Ny() + 2) / mpi_comm_size + 1;
  int64_t Nz_local = 1;

  //int64_t Nx_local = get_Nx();
  //int64_t Ny_local = get_Ny() / mpi_comm_size; //   int64_t Ny_local = (get_Ny() + 2) / mpi_comm_size + 1;
  //int64_t Nz_local = 1;

  //Nx_local_size = get_Nx();
  //Ny_local_size = get_Ny() / mpi_comm_size;

  std::cout << "Ny_Local = " << Ny_local << std::endl;
  Field::ComputingField Ex_local(Nx_local - 2, Ny_local - 2);
  Field::ComputingField Ey_local(Nx_local - 2, Ny_local - 2);
  Field::ComputingField Ez_local(Nx_local - 2, Ny_local - 2);
  Field::ComputingField Bx_local(Nx_local - 2, Ny_local - 2);
  Field::ComputingField By_local(Nx_local - 2, Ny_local - 2);
  Field::ComputingField Bz_local(Nx_local - 2, Ny_local - 2);


  // Пока что предполагается, что Ny кратно количеству процессов size! Также пока что рассматривается MPI для двумерного случая!
  //MPI_Scatter(Ex.data(), Nx_local * Ny_local, MPI_DOUBLE, Ex_local.data(), Nx_local * Ny_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //MPI_Scatter(Ey.data(), Nx_local * Ny_local, MPI_DOUBLE, Ey_local.data(), Nx_local * Ny_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //MPI_Scatter(Ez.data(), Nx_local * Ny_local, MPI_DOUBLE, Ez_local.data(), Nx_local * Ny_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  //MPI_Scatter(Bx.data(), Nx_local * Ny_local, MPI_DOUBLE, Bx_local.data(), Nx_local * Ny_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //MPI_Scatter(By.data(), Nx_local * Ny_local, MPI_DOUBLE, By_local.data(), Nx_local * Ny_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //MPI_Scatter(Bz.data(), Nx_local * Ny_local, MPI_DOUBLE, Bz_local.data(), Nx_local * Ny_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);


  MPI_Scatter(Ex.data(), Ex_local.size(), MPI_DOUBLE, Ex_local.data(), Ex_local.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(Ey.data(), Ey_local.size(), MPI_DOUBLE, Ey_local.data(), Ey_local.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(Ez.data(), Ez_local.size(), MPI_DOUBLE, Ez_local.data(), Ez_local.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Scatter(Bx.data(), Bx_local.size(), MPI_DOUBLE, Bx_local.data(), Bx_local.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(By.data(), By_local.size(), MPI_DOUBLE, By_local.data(), By_local.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(Bz.data(), Bz_local.size(), MPI_DOUBLE, Bz_local.data(), Bz_local.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for (int64_t time = 0ull; time < t; time++)
  {
    //for (i = 0; i < Nx_local; ++i)
    //  for (j = 0; j < Ny_local; ++j)
    //    for (k = 0; k < Nz_local; ++k)
    for (i = 0; i < Nx_local - 2; ++i)
      for (j = 0; j < Ny_local - 2; ++j)
        for (k = 0; k < Nz_local; ++k)
        {
          Bx_local(i, j, k) = Bx_local(i, j, k) + C * B_dt * (/*(Ey_local(i, j, k + 1) - Ey_local(i, j, k)) / dz*/ -(Ez_local(i, j + 1, k) - Ez_local(i, j, k)) / dy);
          By_local(i, j, k) = By_local(i, j, k) + C * B_dt * ((Ez_local(i + 1, j, k) - Ez_local(i, j, k)) / dx /*- (Ex_local(i, j, k + 1) - Ex_local(i, j, k)) / dz*/);
          Bz_local(i, j, k) = Bz_local(i, j, k) + C * B_dt * ((Ex_local(i, j + 1, k) - Ex_local(i, j, k)) / dy - (Ey_local(i + 1, j, k) - Ey_local(i, j, k)) / dx);

        }
    // Обновляю границу - верх 2-х мерной системы
    if (rank == 0 && mpi_comm_size != 1)
    {
      // Получаю верхнюю строчку последнего процесса для синхронизации нулевой строчки нулевого процесса
      // Замечание: здесь get_Nx() = Nx_local - 2
      MPI_Recv(Bx_local.data() + 1, Bx_local.get_Nx(), MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // снизу 
      MPI_Recv(By_local.data() + 1, By_local.get_Nx(), MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // снизу
      MPI_Recv(Bz_local.data() + 1, Bz_local.get_Nx(), MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // снизу

      MPI_Recv(Bx_local.data(), 1, MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // левый нижний узел
      MPI_Recv(By_local.data(), 1, MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // левый нижний узел
      MPI_Recv(Bz_local.data(), 1, MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // левый нижний узел

      MPI_Recv(Bx_local.data() + Bx_local.get_Nx() + 1, 1, MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // правый нижний узел
      MPI_Recv(By_local.data() + By_local.get_Nx() + 1, 1, MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // правый нижний узел
      MPI_Recv(Bz_local.data() + Bz_local.get_Nx() + 1, 1, MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // правый нижний узел


      MPI_Send(&Bx_local(0, 0), Bx_local.get_Nx(), MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // Верх последнего процесса
      MPI_Send(&By_local(0, 0), By_local.get_Nx(), MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // Верх последнего процесса
      MPI_Send(&Bz_local(0, 0), Bz_local.get_Nx(), MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // Верх последнего процесса

      MPI_Send(&Bx_local(0, 0), 1, MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // левый нижний угол (для левого врехнего угла посл. процесса)
      MPI_Send(&By_local(0, 0), 1, MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // левый нижний угол (для левого врехнего угла посл. процесса)
      MPI_Send(&Bz_local(0, 0), 1, MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // левый нижний угол (для левого врехнего угла посл. процесса)

      MPI_Send(&Bx_local(Bx_local.get_Nx() - 1, -1), 1, MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // правый нижний угол (для прав. верх. угла посл. процесса)
      MPI_Send(&By_local(By_local.get_Nx() - 1, -1), 1, MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // правый нижний угол (для прав. верх. угла посл. процесса)
      MPI_Send(&Bz_local(Bz_local.get_Nx() - 1, -1), 1, MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // правый нижний угол (для прав. верх. угла посл. процесса)
    }
    if (rank == mpi_comm_size - 1 && mpi_comm_size != 1)
    {
      // Отправляю верхнюю строчку последнего процесса нулевому процессу для нулевой строчки

      MPI_Send(Bx_local.data() + Nx_local * (Ny_local - 2) + 1, Bx_local.get_Nx(), MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // Низ нулевого процесса
      MPI_Send(By_local.data() + Nx_local * (Ny_local - 2) + 1, By_local.get_Nx(), MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // Низ нулевого процесса
      MPI_Send(Bz_local.data() + Nx_local * (Ny_local - 2) + 1, Bz_local.get_Nx(), MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // Низ нулевого процесса

      MPI_Send(Bx_local.data() + Nx_local * (Ny_local - 1) - 2, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // левый нижний узел нулевого процесса
      MPI_Send(By_local.data() + Nx_local * (Ny_local - 1) - 2, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // левый нижний узел нулевого процесса
      MPI_Send(Bz_local.data() + Nx_local * (Ny_local - 1) - 2, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // левый нижний узел нулевого процесса

      MPI_Send(Bx_local.data() + Nx_local * (Ny_local - 2) + 1, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // правый нижний узел узел нулевого процесса
      MPI_Send(By_local.data() + Nx_local * (Ny_local - 2) + 1, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // правый нижний узел узел нулевого процесса
      MPI_Send(Bz_local.data() + Nx_local * (Ny_local - 2) + 1, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // правый нижний узел узел нулевого процесса

      MPI_Recv(&Bx_local(0, Bx_local.get_Ny() - 2), Bx_local.get_Nx(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // сверху (получение нулевой строки нулевого процесса)
      MPI_Recv(&By_local(0, By_local.get_Ny() - 2), By_local.get_Nx(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // сверху (получение нулевой строки нулевого процесса)
      MPI_Recv(&Bz_local(0, Bz_local.get_Ny() - 2), Bz_local.get_Nx(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // сверху (получение нулевой строки нулевого процесса)

      MPI_Recv(&Bx_local(Bx_local.get_Nx(), Bx_local.get_Ny()), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // правый верхний узел (получение из нулевого процесса)
      MPI_Recv(&By_local(By_local.get_Nx(), By_local.get_Ny()), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // правый верхний узел (получение из нулевого процесса)
      MPI_Recv(&Bz_local(Bz_local.get_Nx(), Bz_local.get_Ny()), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // правый верхний узел (получение из нулевого процесса)

      MPI_Recv(&Bx_local(-1, Bx_local.get_Ny() - 1), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // левый верхний узел (получение из нулевого процесса)
      MPI_Recv(&By_local(-1, By_local.get_Ny() - 1), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // левый верхний узел (получение из нулевого процесса)
      MPI_Recv(&Bz_local(-1, Bz_local.get_Ny() - 1), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // левый верхний узел (получение из нулевого процесса)

    }


    if (mpi_comm_size == 1)
    {
      for (int64_t x = 0; x < Bx_local.get_Nx(); ++x)
      {
        *(Bx_local.data() + x + 1) = Bx_local(x, Bx_local.get_Ny() - 1); // снизу
        Bx_local(x, Bx_local.get_Ny()) = Bx_local(x, 0); // сверху

        *(By_local.data() + x + 1) = By_local(x, By_local.get_Ny() - 1); // снизу
        By_local(x, By_local.get_Ny()) = By_local(x, 0); // сверху

        *(Bz_local.data() + x + 1) = Bz_local(x, Bz_local.get_Ny() - 1); // снизу
        Bz_local(x, Bz_local.get_Ny()) = Bz_local(x, 0); // сверху
      }
      for (int64_t y = 0; y < Bx_local.get_Ny(); ++y)
      {
        Bx_local(-1, y) = Bx_local(Bx_local.get_Nx() - 1, y); // слева
        Bx_local(Bx_local.get_Nx(), y) = Bx_local(0, y); // справа

        By_local(-1, y) = By_local(By_local.get_Nx() - 1, y); // слева
        By_local(By_local.get_Nx(), y) = By_local(0, y); // справа

        Bz_local(-1, y) = Bz_local(Bz_local.get_Nx() - 1, y); // слева
        Bz_local(Bz_local.get_Nx(), y) = Bz_local(0, y); // справа
      }
      Bx_local(-1, -1) = Bx_local(Bx_local.get_Nx() - 1, Bx_local.get_Ny() - 1); // левый нижний узел
      Bx_local(-1, Bx_local.get_Ny()) = Bx_local(Bx_local.get_Nx() - 1, 0); // левый верхний узел
      Bx_local(Bx_local.get_Nx(), Bx_local.get_Ny()) = Bx_local(0, 0); // правый верхний узел
      *(Bx_local.data() + Bx_local.get_Nx() + 1) = Bx_local(0, Bx_local.get_Ny() - 1); // правый нижний узел

      By_local(-1, -1) = By_local(By_local.get_Nx() - 1, By_local.get_Ny() - 1); // левый нижний узел
      By_local(-1, By_local.get_Ny()) = By_local(By_local.get_Nx() - 1, 0); // левый верхний узел
      By_local(By_local.get_Nx(), By_local.get_Ny()) = By_local(0, 0); // правый верхний узел
      *(By_local.data() + By_local.get_Nx() + 1) = By_local(0, By_local.get_Ny() - 1); // правый нижний узел

      Bz_local(-1, -1) = Bz_local(Bz_local.get_Nx() - 1, Bz_local.get_Ny() - 1); // левый нижний узел
      Bz_local(-1, Bz_local.get_Ny()) = Bz_local(Bz_local.get_Nx() - 1, 0); // левый верхний узел
      Bz_local(Bz_local.get_Nx(), Bz_local.get_Ny()) = Bz_local(0, 0); // правый верхний узел
      *(Bz_local.data() + Bz_local.get_Nx() + 1) = Bz_local(0, Bz_local.get_Ny() - 1); // правый нижний узел
    }


    // Синхронизация смежных строк (для 2-мерного случая)
    if (mpi_comm_size > 1)
    {
      if (rank != mpi_comm_size - 1 && rank != 0)
      {
        MPI_Sendrecv(Bx_local.data() + Nx_local * (Ny_local - 2), Nx_local, MPI_DOUBLE, rank + 1, rank,
          Bx_local.data() + Nx_local * (Ny_local - 1), Nx_local, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status); // сверху получил-отправил
        MPI_Sendrecv(By_local.data() + Nx_local * (Ny_local - 2), Nx_local, MPI_DOUBLE, rank + 1, rank,
          By_local.data() + Nx_local * (Ny_local - 1), Nx_local, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status); // сверху получил-отправил
        MPI_Sendrecv(Bz_local.data() + Nx_local * (Ny_local - 2), Nx_local, MPI_DOUBLE, rank + 1, rank,
          Bz_local.data() + Nx_local * (Ny_local - 1), Nx_local, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status); // сверху получил-отправил

        MPI_Sendrecv(Bx_local.data() + Nx_local, Nx_local, MPI_DOUBLE, rank - 1, rank,
          Bx_local.data(), Nx_local, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status); // снизу получил-отправил
        MPI_Sendrecv(By_local.data() + Nx_local, Nx_local, MPI_DOUBLE, rank - 1, rank,
          By_local.data(), Nx_local, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status); // снизу получил-отправил
        MPI_Sendrecv(Bz_local.data() + Nx_local, Nx_local, MPI_DOUBLE, rank - 1, rank,
          Bz_local.data(), Nx_local, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status); // снизу получил-отправил
      }
      else if (rank == 0)
      {
        MPI_Sendrecv(Bx_local.data() + Nx_local * (Ny_local - 2), Nx_local, MPI_DOUBLE, rank + 1, rank,
          Bx_local.data() + Nx_local * (Ny_local - 1), Nx_local, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status); // сверху получил-отправил
        MPI_Sendrecv(By_local.data() + Nx_local * (Ny_local - 2), Nx_local, MPI_DOUBLE, rank + 1, rank,
          By_local.data() + Nx_local * (Ny_local - 1), Nx_local, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status); // сверху получил-отправил
        MPI_Sendrecv(Bz_local.data() + Nx_local * (Ny_local - 2), Nx_local, MPI_DOUBLE, rank + 1, rank,
          Bz_local.data() + Nx_local * (Ny_local - 1), Nx_local, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status); // сверху получил-отправил
      }
      else // if (rank == mpi_comm_size - 1)
      {
        MPI_Sendrecv(Bx_local.data() + Nx_local, Nx_local, MPI_DOUBLE, rank - 1, rank,
          Bx_local.data(), Nx_local, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status); // снизу получил-отправил
        MPI_Sendrecv(By_local.data() + Nx_local, Nx_local, MPI_DOUBLE, rank - 1, rank,
          By_local.data(), Nx_local, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status); // снизу получил-отправил
        MPI_Sendrecv(Bz_local.data() + Nx_local, Nx_local, MPI_DOUBLE, rank - 1, rank,
          Bz_local.data(), Nx_local, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status); // снизу получил-отправил
      }
    }
// ==============================================================================================================================================

    //for (i = 0ull; i < Nx_local; ++i)
    //  for (j = 0ull; j < Ny_local; ++j)
    //    for (k = 0ull; k < Nz_local; ++k)
    for (i = 0; i < Nx_local - 2; ++i)
      for (j = 0; j < Ny_local - 2; ++j)
        for (k = 0; k < Nz_local; ++k)
        {
          Ex_local(i, j, k) = Ex_local(i, j, k) + C * E_dt * ((Bz_local(i, j, k) - Bz_local(i, j - 1, k)) / dy /*- (By_local(i, j, k) - By_local(i, j, k - 1)) / dz*/);
          Ey_local(i, j, k) = Ey_local(i, j, k) + C * E_dt * (/*(Bx_local(i, j, k) - Bx_local(i, j, k - 1)) / dz*/ - (Bz_local(i, j, k) - Bz_local(i - 1, j, k)) / dx);
          Ez_local(i, j, k) = Ez_local(i, j, k) + C * E_dt * ((By_local(i, j, k) - By_local(i - 1, j, k)) / dx - (Bx_local(i, j, k) - Bx_local(i, j - 1, k)) / dy);
        }
    // Обновляю границу - верх 2-х мерной системы
    if (rank == 0 && mpi_comm_size != 1)
    {
      // Получаю верхнюю строчку последнего процесса для синхронизации нулевой строчки нулевого процесса
      // Замечание: здесь get_Nx() = Nx_local - 2
      MPI_Recv(Ex_local.data() + 1, Ex_local.get_Nx(), MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // снизу 
      MPI_Recv(Ey_local.data() + 1, Ey_local.get_Nx(), MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // снизу
      MPI_Recv(Ez_local.data() + 1, Ez_local.get_Nx(), MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // снизу

      MPI_Recv(Ex_local.data(), 1, MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // левый нижний узел
      MPI_Recv(Ey_local.data(), 1, MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // левый нижний узел
      MPI_Recv(Ez_local.data(), 1, MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // левый нижний узел

      MPI_Recv(Ex_local.data() + Ex_local.get_Nx() + 1, 1, MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // правый нижний узел
      MPI_Recv(Ey_local.data() + Ey_local.get_Nx() + 1, 1, MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // правый нижний узел
      MPI_Recv(Ez_local.data() + Ez_local.get_Nx() + 1, 1, MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // правый нижний узел

      MPI_Send(&Ex_local(0, 0), Ex_local.get_Nx(), MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // Верх последнего процесса
      MPI_Send(&Ey_local(0, 0), Ey_local.get_Nx(), MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // Верх последнего процесса
      MPI_Send(&Ez_local(0, 0), Ez_local.get_Nx(), MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // Верх последнего процесса

      MPI_Send(&Ex_local(0, 0), 1, MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // левый нижний угол (для левого врехнего угла посл. процесса)
      MPI_Send(&Ey_local(0, 0), 1, MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // левый нижний угол (для левого врехнего угла посл. процесса)
      MPI_Send(&Ez_local(0, 0), 1, MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // левый нижний угол (для левого врехнего угла посл. процесса)

      MPI_Send(&Ex_local(Ex_local.get_Nx() - 1, -1), 1, MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // правый нижний угол (для прав. верх. угла посл. процесса)
      MPI_Send(&Ey_local(Ey_local.get_Nx() - 1, -1), 1, MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // правый нижний угол (для прав. верх. угла посл. процесса)
      MPI_Send(&Ez_local(Ez_local.get_Nx() - 1, -1), 1, MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // правый нижний угол (для прав. верх. угла посл. процесса)
    }
    if (rank == mpi_comm_size - 1 && mpi_comm_size != 1)
    {
      // Отправляю верхнюю строчку последнего процесса нулевому процессу для нулевой строчки

      MPI_Send(Ex_local.data() + Nx_local * (Ny_local - 2) + 1, Ex_local.get_Nx(), MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // Низ нулевого процесса
      MPI_Send(Ey_local.data() + Nx_local * (Ny_local - 2) + 1, Ey_local.get_Nx(), MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // Низ нулевого процесса
      MPI_Send(Ez_local.data() + Nx_local * (Ny_local - 2) + 1, Ez_local.get_Nx(), MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // Низ нулевого процесса

      MPI_Send(Ex_local.data() + Nx_local * (Ny_local - 1) - 2, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // левый нижний узел нулевого процесса
      MPI_Send(Ey_local.data() + Nx_local * (Ny_local - 1) - 2, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // левый нижний узел нулевого процесса
      MPI_Send(Ez_local.data() + Nx_local * (Ny_local - 1) - 2, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // левый нижний узел нулевого процесса

      MPI_Send(Ex_local.data() + Nx_local * (Ny_local - 2) + 1, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // правый нижний узел узел нулевого процесса
      MPI_Send(Ey_local.data() + Nx_local * (Ny_local - 2) + 1, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // правый нижний узел узел нулевого процесса
      MPI_Send(Ez_local.data() + Nx_local * (Ny_local - 2) + 1, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // правый нижний узел узел нулевого процесса

      MPI_Recv(&Ex_local(0, Ex_local.get_Ny() - 2), Ex_local.get_Nx(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // сверху (получение нулевой строки нулевого процесса)
      MPI_Recv(&Ey_local(0, Ey_local.get_Ny() - 2), Ey_local.get_Nx(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // сверху (получение нулевой строки нулевого процесса)
      MPI_Recv(&Ez_local(0, Ez_local.get_Ny() - 2), Ez_local.get_Nx(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // сверху (получение нулевой строки нулевого процесса)

      MPI_Recv(&Ex_local(Ex_local.get_Nx(), Ex_local.get_Ny()), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // правый верхний узел (получение из нулевого процесса)
      MPI_Recv(&Ey_local(Ey_local.get_Nx(), Ey_local.get_Ny()), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // правый верхний узел (получение из нулевого процесса)
      MPI_Recv(&Ez_local(Ez_local.get_Nx(), Ez_local.get_Ny()), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // правый верхний узел (получение из нулевого процесса)

      MPI_Recv(&Ex_local(-1, Ex_local.get_Ny() - 1), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // левый верхний узел (получение из нулевого процесса)
      MPI_Recv(&Ey_local(-1, Ey_local.get_Ny() - 1), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // левый верхний узел (получение из нулевого процесса)
      MPI_Recv(&Ez_local(-1, Ez_local.get_Ny() - 1), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // левый верхний узел (получение из нулевого процесса)
    }

    if (mpi_comm_size == 1)
    {
      for (int64_t x = 0; x < Ex_local.get_Nx(); ++x)
      {
        *(Ex_local.data() + x + 1) = Ex_local(x, Ex_local.get_Ny() - 1); // снизу
        Ex_local(x, Ex_local.get_Ny()) = Ex_local(x, 0); // сверху
    
        *(Ey_local.data() + x + 1) = Ey_local(x, Ey_local.get_Ny() - 1); // снизу
        Ey_local(x, Ey_local.get_Ny()) = Ey_local(x, 0); // сверху
    
        *(Ez_local.data() + x + 1) = Ez_local(x, Ez_local.get_Ny() - 1); // снизу
        Ez_local(x, Ez_local.get_Ny()) = Ez_local(x, 0); // сверху
      }
      for (int64_t y = 0; y < Ex_local.get_Ny(); ++y)
      {
        Ex_local(-1, y) = Ex_local(Ex_local.get_Nx() - 1, y); // слева
        Ex_local(Ex_local.get_Nx(), y) = Ex_local(0, y); // справа
    
        Ey_local(-1, y) = Ey_local(Ey_local.get_Nx() - 1, y); // слева
        Ey_local(Ey_local.get_Nx(), y) = Ey_local(0, y); // справа
    
        Ez_local(-1, y) = Ez_local(Ez_local.get_Nx() - 1, y); // слева
        Ez_local(Ez_local.get_Nx(), y) = Ez_local(0, y); // справа
      }
      Ex_local(-1, -1) = Ex_local(Ex_local.get_Nx() - 1, Ex_local.get_Ny() - 1); // левый нижний узел
      Ex_local(-1, Ex_local.get_Ny()) = Ex_local(Ex_local.get_Nx() - 1, 0); // левый верхний узел
      Ex_local(Ex_local.get_Nx(), Ex_local.get_Ny()) = Ex_local(0, 0); // правый верхний узел
      *(Ex_local.data() + Ex_local.get_Nx() + 1) = Ex_local(0, Ex_local.get_Ny() - 1); // правый нижний узел
    
      Ey_local(-1, -1) = Ey_local(Ey_local.get_Nx() - 1, Ey_local.get_Ny() - 1); // левый нижний узел
      Ey_local(-1, Ey_local.get_Ny()) = Ey_local(Ey_local.get_Nx() - 1, 0); // левый верхний узел
      Ey_local(Ey_local.get_Nx(), Ey_local.get_Ny()) = Ey_local(0, 0); // правый верхний узел
      *(Ey_local.data() + Ey_local.get_Nx() + 1) = Ey_local(0, Ey_local.get_Ny() - 1); // правый нижний узел
    
      Ez_local(-1, -1) = Ez_local(Ez_local.get_Nx() - 1, Ez_local.get_Ny() - 1); // левый нижний узел
      Ez_local(-1, Ez_local.get_Ny()) = Ez_local(Ez_local.get_Nx() - 1, 0); // левый верхний узел
      Ez_local(Ez_local.get_Nx(), Ez_local.get_Ny()) = Ez_local(0, 0); // правый верхний узел
      *(Ez_local.data() + Ez_local.get_Nx() + 1) = Ez_local(0, Ez_local.get_Ny() - 1); // правый нижний узел
    }

    // Синхронизация смежных строк (для 2-мерного случая)
    if (mpi_comm_size > 1)
    {
      if (rank != mpi_comm_size - 1 && rank != 0)
      {
        MPI_Sendrecv(Ex_local.data() + Nx_local * (Ny_local - 2), Nx_local, MPI_DOUBLE, rank + 1, rank,
          Ex_local.data() + Nx_local * (Ny_local - 1), Nx_local, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status); // сверху получил-отправил
        MPI_Sendrecv(Ey_local.data() + Nx_local * (Ny_local - 2), Nx_local, MPI_DOUBLE, rank + 1, rank,
          Ey_local.data() + Nx_local * (Ny_local - 1), Nx_local, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status); // сверху получил-отправил
        MPI_Sendrecv(Ez_local.data() + Nx_local * (Ny_local - 2), Nx_local, MPI_DOUBLE, rank + 1, rank,
          Ez_local.data() + Nx_local * (Ny_local - 1), Nx_local, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status); // сверху получил-отправил

        MPI_Sendrecv(Ex_local.data() + Nx_local, Nx_local, MPI_DOUBLE, rank - 1, rank,
          Ex_local.data(), Nx_local, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status); // снизу получил-отправил
        MPI_Sendrecv(Ey_local.data() + Nx_local, Nx_local, MPI_DOUBLE, rank - 1, rank,
          Ey_local.data(), Nx_local, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status); // снизу получил-отправил
        MPI_Sendrecv(Ez_local.data() + Nx_local, Nx_local, MPI_DOUBLE, rank - 1, rank,
          Ez_local.data(), Nx_local, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status); // снизу получил-отправил
      }
      else if (rank == 0)
      {
        MPI_Sendrecv(Ex_local.data() + Nx_local * (Ny_local - 2), Nx_local, MPI_DOUBLE, rank + 1, rank,
          Ex_local.data() + Nx_local * (Ny_local - 1), Nx_local, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status); // сверху получил-отправил
        MPI_Sendrecv(Ey_local.data() + Nx_local * (Ny_local - 2), Nx_local, MPI_DOUBLE, rank + 1, rank,
          Ey_local.data() + Nx_local * (Ny_local - 1), Nx_local, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status); // сверху получил-отправил
        MPI_Sendrecv(Ez_local.data() + Nx_local * (Ny_local - 2), Nx_local, MPI_DOUBLE, rank + 1, rank,
          Ez_local.data() + Nx_local * (Ny_local - 1), Nx_local, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status); // сверху получил-отправил
      }
      else // if (rank == mpi_comm_size - 1)
      {
        MPI_Sendrecv(Ex_local.data() + Nx_local, Nx_local, MPI_DOUBLE, rank - 1, rank,
          Ex_local.data(), Nx_local, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status); // снизу получил-отправил
        MPI_Sendrecv(Ey_local.data() + Nx_local, Nx_local, MPI_DOUBLE, rank - 1, rank,
          Ey_local.data(), Nx_local, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status); // снизу получил-отправил
        MPI_Sendrecv(Ez_local.data() + Nx_local, Nx_local, MPI_DOUBLE, rank - 1, rank,
          Ez_local.data(), Nx_local, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status); // снизу получил-отправил
      }
    }

// ==============================================================================================================================================

    //for (i = 0; i < Nx_local; ++i)
    //  for (j = 0; j < Ny_local; ++j)
    //    for (k = 0; k < Nz_local; ++k)

    for (i = 0; i < Nx_local - 2; ++i)
      for (j = 0; j < Ny_local - 2; ++j)
        for (k = 0; k < Nz_local; ++k)
        {
          Bx_local(i, j, k) = Bx_local(i, j, k) + C * B_dt * (/*(Ey_local(i, j, k + 1) - Ey_local(i, j, k)) / dz*/ - (Ez_local(i, j + 1, k) - Ez_local(i, j, k)) / dy);
          By_local(i, j, k) = By_local(i, j, k) + C * B_dt * ((Ez_local(i + 1, j, k) - Ez_local(i, j, k)) / dx /*- (Ex_local(i, j, k + 1) - Ex_local(i, j, k)) / dz*/);
          Bz_local(i, j, k) = Bz_local(i, j, k) + C * B_dt * ((Ex_local(i, j + 1, k) - Ex_local(i, j, k)) / dy - (Ey_local(i + 1, j, k) - Ey_local(i, j, k)) / dx);
          //std::cout << "i = " << i << " j = " << j << " k = " << k << '\n';
        }
    // Обновляю границу - верх 2-х мерной системы
    if (rank == 0 && mpi_comm_size != 1)
    {
      // Получаю верхнюю строчку последнего процесса для синхронизации нулевой строчки нулевого процесса
      // Замечание: здесь get_Nx() = Nx_local - 2
      MPI_Recv(Bx_local.data() + 1, Bx_local.get_Nx(), MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // снизу 
      MPI_Recv(By_local.data() + 1, By_local.get_Nx(), MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // снизу
      MPI_Recv(Bz_local.data() + 1, Bz_local.get_Nx(), MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // снизу

      MPI_Recv(Bx_local.data(), 1, MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // левый нижний узел
      MPI_Recv(By_local.data(), 1, MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // левый нижний узел
      MPI_Recv(Bz_local.data(), 1, MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // левый нижний узел

      MPI_Recv(Bx_local.data() + Bx_local.get_Nx() + 1, 1, MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // правый нижний узел
      MPI_Recv(By_local.data() + By_local.get_Nx() + 1, 1, MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // правый нижний узел
      MPI_Recv(Bz_local.data() + Bz_local.get_Nx() + 1, 1, MPI_DOUBLE, mpi_comm_size - 1, mpi_comm_size - 1, MPI_COMM_WORLD, &status); // правый нижний узел

      MPI_Send(&Bx_local(0, 0), Bx_local.get_Nx(), MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // Верх последнего процесса
      MPI_Send(&By_local(0, 0), By_local.get_Nx(), MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // Верх последнего процесса
      MPI_Send(&Bz_local(0, 0), Bz_local.get_Nx(), MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // Верх последнего процесса

      MPI_Send(&Bx_local(0, 0), 1, MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // левый нижний угол (для левого врехнего угла посл. процесса)
      MPI_Send(&By_local(0, 0), 1, MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // левый нижний угол (для левого врехнего угла посл. процесса)
      MPI_Send(&Bz_local(0, 0), 1, MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // левый нижний угол (для левого врехнего угла посл. процесса)

      MPI_Send(&Bx_local(Bx_local.get_Nx() - 1, -1), 1, MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // правый нижний угол (для прав. верх. угла посл. процесса)
      MPI_Send(&By_local(By_local.get_Nx() - 1, -1), 1, MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // правый нижний угол (для прав. верх. угла посл. процесса)
      MPI_Send(&Bz_local(Bz_local.get_Nx() - 1, -1), 1, MPI_DOUBLE, mpi_comm_size - 1, rank, MPI_COMM_WORLD); // правый нижний угол (для прав. верх. угла посл. процесса)
    }
    if (rank == mpi_comm_size - 1 && mpi_comm_size != 1)
    {
      // Отправляю верхнюю строчку последнего процесса нулевому процессу для нулевой строчки

      MPI_Send(Bx_local.data() + Nx_local * (Ny_local - 2) + 1, Bx_local.get_Nx(), MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // Низ нулевого процесса
      MPI_Send(By_local.data() + Nx_local * (Ny_local - 2) + 1, By_local.get_Nx(), MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // Низ нулевого процесса
      MPI_Send(Bz_local.data() + Nx_local * (Ny_local - 2) + 1, Bz_local.get_Nx(), MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // Низ нулевого процесса

      MPI_Send(Bx_local.data() + Nx_local * (Ny_local - 1) - 2, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // левый нижний узел нулевого процесса
      MPI_Send(By_local.data() + Nx_local * (Ny_local - 1) - 2, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // левый нижний узел нулевого процесса
      MPI_Send(Bz_local.data() + Nx_local * (Ny_local - 1) - 2, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // левый нижний узел нулевого процесса

      MPI_Send(Bx_local.data() + Nx_local * (Ny_local - 2) + 1, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // правый нижний узел узел нулевого процесса
      MPI_Send(By_local.data() + Nx_local * (Ny_local - 2) + 1, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // правый нижний узел узел нулевого процесса
      MPI_Send(Bz_local.data() + Nx_local * (Ny_local - 2) + 1, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); // правый нижний узел узел нулевого процесса

      MPI_Recv(&Bx_local(0, Bx_local.get_Ny() - 2), Bx_local.get_Nx(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // сверху (получение нулевой строки нулевого процесса)
      MPI_Recv(&By_local(0, By_local.get_Ny() - 2), By_local.get_Nx(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // сверху (получение нулевой строки нулевого процесса)
      MPI_Recv(&Bz_local(0, Bz_local.get_Ny() - 2), Bz_local.get_Nx(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // сверху (получение нулевой строки нулевого процесса)

      MPI_Recv(&Bx_local(Bx_local.get_Nx(), Bx_local.get_Ny()), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // правый верхний узел (получение из нулевого процесса)
      MPI_Recv(&By_local(By_local.get_Nx(), By_local.get_Ny()), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // правый верхний узел (получение из нулевого процесса)
      MPI_Recv(&Bz_local(Bz_local.get_Nx(), Bz_local.get_Ny()), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // правый верхний узел (получение из нулевого процесса)

      MPI_Recv(&Bx_local(-1, Bx_local.get_Ny() - 1), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // левый верхний узел (получение из нулевого процесса)
      MPI_Recv(&By_local(-1, By_local.get_Ny() - 1), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // левый верхний узел (получение из нулевого процесса)
      MPI_Recv(&Bz_local(-1, Bz_local.get_Ny() - 1), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // левый верхний узел (получение из нулевого процесса)
    }

    // Синхронизация смежных строк (для 2-мерного случая)
    if (mpi_comm_size > 1)
    {
      if (rank != mpi_comm_size - 1 && rank != 0)
      {
        MPI_Sendrecv(Bx_local.data() + Nx_local * (Ny_local - 2), Nx_local, MPI_DOUBLE, rank + 1, rank,
          Bx_local.data() + Nx_local * (Ny_local - 1), Nx_local, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status); // сверху получил-отправил
        MPI_Sendrecv(By_local.data() + Nx_local * (Ny_local - 2), Nx_local, MPI_DOUBLE, rank + 1, rank,
          By_local.data() + Nx_local * (Ny_local - 1), Nx_local, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status); // сверху получил-отправил
        MPI_Sendrecv(Bz_local.data() + Nx_local * (Ny_local - 2), Nx_local, MPI_DOUBLE, rank + 1, rank,
          Bz_local.data() + Nx_local * (Ny_local - 1), Nx_local, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status); // сверху получил-отправил

        MPI_Sendrecv(Bx_local.data() + Nx_local, Nx_local, MPI_DOUBLE, rank - 1, rank,
          Bx_local.data(), Nx_local, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status); // снизу получил-отправил
        MPI_Sendrecv(By_local.data() + Nx_local, Nx_local, MPI_DOUBLE, rank - 1, rank,
          By_local.data(), Nx_local, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status); // снизу получил-отправил
        MPI_Sendrecv(Bz_local.data() + Nx_local, Nx_local, MPI_DOUBLE, rank - 1, rank,
          Bz_local.data(), Nx_local, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status); // снизу получил-отправил
      }
      else if (rank == 0)
      {
        MPI_Sendrecv(Bx_local.data() + Nx_local * (Ny_local - 2), Nx_local, MPI_DOUBLE, rank + 1, rank,
          Bx_local.data() + Nx_local * (Ny_local - 1), Nx_local, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status); // сверху получил-отправил
        MPI_Sendrecv(By_local.data() + Nx_local * (Ny_local - 2), Nx_local, MPI_DOUBLE, rank + 1, rank,
          By_local.data() + Nx_local * (Ny_local - 1), Nx_local, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status); // сверху получил-отправил
        MPI_Sendrecv(Bz_local.data() + Nx_local * (Ny_local - 2), Nx_local, MPI_DOUBLE, rank + 1, rank,
          Bz_local.data() + Nx_local * (Ny_local - 1), Nx_local, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status); // сверху получил-отправил
      }
      else // if (rank == mpi_comm_size - 1)
      {
        MPI_Sendrecv(Bx_local.data() + Nx_local, Nx_local, MPI_DOUBLE, rank - 1, rank,
          Bx_local.data(), Nx_local, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status); // снизу получил-отправил
        MPI_Sendrecv(By_local.data() + Nx_local, Nx_local, MPI_DOUBLE, rank - 1, rank,
          By_local.data(), Nx_local, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status); // снизу получил-отправил
        MPI_Sendrecv(Bz_local.data() + Nx_local, Nx_local, MPI_DOUBLE, rank - 1, rank,
          Bz_local.data(), Nx_local, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status); // снизу получил-отправил
      }
    }


  }

  std::vector<double> proc_buffer_Ex(Nx_local* (Ny_local - 1));
  std::vector<double> proc_buffer_Ey(Nx_local* (Ny_local - 1));
  std::vector<double> proc_buffer_Ez(Nx_local* (Ny_local - 1));

  std::vector<double> proc_buffer_Bx(Nx_local* (Ny_local - 1));
  std::vector<double> proc_buffer_By(Nx_local* (Ny_local - 1));
  std::vector<double> proc_buffer_Bz(Nx_local* (Ny_local - 1));


  std::cout << "proc_buf_Ex size =  " << proc_buffer_Ex.size() << std::endl;
  if (rank == 0)
  {
    for (int i = 0; i < Bx_local.size(); ++i)
    {
      std::cout << *(Bx_local.data() + i) << std::endl;
    }
  }

  // Gather для нулевого процесса
  if (rank == 0)
  {
    memcpy(Ex.data(), Ex_local.data(), proc_buffer_Ex.size() * sizeof(double));
    memcpy(Ey.data(), Ey_local.data(), proc_buffer_Ey.size() * sizeof(double));
    memcpy(Ez.data(), Ez_local.data(), proc_buffer_Ez.size() * sizeof(double));

    memcpy(Bx.data(), Bx_local.data(), proc_buffer_Bx.size() * sizeof(double));
    memcpy(By.data(), By_local.data(), proc_buffer_By.size() * sizeof(double));
    memcpy(Bz.data(), Bz_local.data(), proc_buffer_Bz.size() * sizeof(double));
  }
  else
  {
    memcpy(proc_buffer_Ex.data(), Ex_local.data() + Ex_local.get_Nx(), proc_buffer_Ex.size() * sizeof(double));
    memcpy(proc_buffer_Ey.data(), Ey_local.data() + Ey_local.get_Nx(), proc_buffer_Ey.size() * sizeof(double));
    memcpy(proc_buffer_Ez.data(), Ez_local.data() + Ez_local.get_Nx(), proc_buffer_Ez.size() * sizeof(double));
    
    memcpy(proc_buffer_Bx.data(), Bx_local.data() + Bx_local.get_Nx(), proc_buffer_Bx.size() * sizeof(double));
    memcpy(proc_buffer_By.data(), By_local.data() + By_local.get_Nx(), proc_buffer_By.size() * sizeof(double));
    memcpy(proc_buffer_Bz.data(), Bz_local.data() + Bz_local.get_Nx(), proc_buffer_Bz.size() * sizeof(double));

    MPI_Gather(proc_buffer_Ex.data(), proc_buffer_Ex.size(), MPI_DOUBLE, Ex.data(), proc_buffer_Ex.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(proc_buffer_Ey.data(), proc_buffer_Ey.size(), MPI_DOUBLE, Ey.data(), proc_buffer_Ey.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(proc_buffer_Ez.data(), proc_buffer_Ez.size(), MPI_DOUBLE, Ez.data(), proc_buffer_Ez.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Gather(proc_buffer_Bx.data(), proc_buffer_Bx.size(), MPI_DOUBLE, Bx.data(), proc_buffer_Bx.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(proc_buffer_By.data(), proc_buffer_By.size(), MPI_DOUBLE, By.data(), proc_buffer_By.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(proc_buffer_Bz.data(), proc_buffer_Bz.size(), MPI_DOUBLE, Bz.data(), proc_buffer_Bz.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  }

    //MPI_Gather(Ex_local.data(), Ex.size(), MPI_DOUBLE, Ex.data(), Ex_local.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Gather(Ey_local.data(), Ey.size(), MPI_DOUBLE, Ey.data(), Ey_local.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Gather(Ez_local.data(), Ez.size(), MPI_DOUBLE, Ez.data(), Ez_local.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
   
    //MPI_Gather(Bx_local.data(), Bx.size(), MPI_DOUBLE, Bx.data(), Bx_local.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Gather(By_local.data(), By.size(), MPI_DOUBLE, By.data(), By_local.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Gather(Bz_local.data(), Bz.size(), MPI_DOUBLE, Bz.data(), Bz_local.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (mpi_comm_size == 1)
    {
      std::memcpy(Ex.data(), Ex_local.data(), Ex_local.size() * sizeof(double));
      std::memcpy(Ey.data(), Ey_local.data(), Ey_local.size() * sizeof(double));
      std::memcpy(Ez.data(), Ez_local.data(), Ez_local.size() * sizeof(double));

      std::memcpy(Bx.data(), Bx_local.data(), Bx_local.size() * sizeof(double));
      std::memcpy(By.data(), By_local.data(), By_local.size() * sizeof(double));
      std::memcpy(Bz.data(), Bz_local.data(), Bz_local.size() * sizeof(double));

      //MPI_Send(Ex_local.data(), Nx_local * Ny_local, MPI_DOUBLE)
    }

//#else  
  //for (int64_t time = 0ull; time < t; time++)
  //{
  //  for (i = 0; i < Nx; ++i)
  //    for (j = 0; j < Ny; ++j)
  //      for (k = 0; k < Nz; ++k)
  //      {
  //        Bx(i, j, k) = Bx(i, j, k) + C * B_dt * (/*(Ey(i, j, k + 1) - Ey(i, j, k)) / dz*/ - (Ez(i, j + 1, k) - Ez(i, j, k)) / dy);
  //        By(i, j, k) = By(i, j, k) + C * B_dt * ((Ez(i + 1, j, k) - Ez(i, j, k)) / dx /*- (Ex(i, j, k + 1) - Ex(i, j, k)) / dz*/);
  //        Bz(i, j, k) = Bz(i, j, k) + C * B_dt * ((Ex(i, j + 1, k) - Ex(i, j, k)) / dy - (Ey(i + 1, j, k) - Ey(i, j, k)) / dx);
  //        //std::cout << "i = " << i << " j = " << j << " k = " << k << '\n';
  //      }
  //  // Обновляю границу - верх 2-х мерной системы
  //  for (int64_t x = 0; x < Bx.get_Nx(); ++x)
  //  {
  //    *(Bx.data() + x + 1) = Bx(x, Bx.get_Ny() - 1); // снизу
  //    Bx(x, Bx.get_Ny()) = Bx(x, 0); // сверху

  //    *(By.data() + x + 1) = By(x, By.get_Ny() - 1); // снизу
  //    By(x, By.get_Ny() - 1) = By(x, 0); // сверху

  //    *(Bz.data() + x + 1) = Bz(x, Bz.get_Ny() - 1); // снизу
  //    Bz(x, Bz.get_Ny()) = Bz(x, 0); // сверху
  //  }
  //  for (int64_t y = 0; y < Bx.get_Ny(); ++y)
  //  {
  //    Bx(-1, y) = Bx(Bx.get_Nx() - 1, y); // слева
  //    Bx(Bx.get_Nx(), y) = Bx(0, y); // справа

  //    By(-1, y) = By(By.get_Nx() - 1, y); // слева
  //    By(By.get_Nx(), y) = By(0, y); // справа

  //    Bz(-1, y) = Bz(Bz.get_Nx() - 1, y); // слева
  //    Bz(Bz.get_Nx(), y) = Bz(0, y); // справа
  //  }
  //  Bx(-1, -1) = Bx(Bx.get_Nx() - 1, Bx.get_Ny() - 1); // левый нижний узел
  //  Bx(-1, Bx.get_Ny()) = Bx(Bx.get_Nx() - 1, 0); // левый верхний узел
  //  Bx(Bx.get_Nx(), Bx.get_Ny()) = Bx(0, 0); // правый верхний узел
  //  *(Bx.data() + Bx.get_Nx() + 1) = Bx(0, Bx.get_Ny() - 1); // правый нижний узел

  //  By(-1, -1) = By(By.get_Nx() - 1, By.get_Ny() - 1); // левый нижний узел
  //  By(-1, By.get_Ny()) = By(By.get_Nx() - 1, 0); // левый верхний узел
  //  By(By.get_Nx(), By.get_Ny()) = By(0, 0); // правый верхний узел
  //  *(By.data() + By.get_Nx() + 1) = By(0, By.get_Ny() - 1); // правый нижний узел

  //  Bz(-1, -1) = Bz(Bz.get_Nx() - 1, Bz.get_Ny() - 1); // левый нижний узел
  //  Bz(-1, Bz.get_Ny()) = Bz(Bz.get_Nx() - 1, 0); // левый верхний узел
  //  Bz(Bz.get_Nx(), Bz.get_Ny()) = Bz(0, 0); // правый верхний узел
  //  *(Bz.data() + Bz.get_Nx() + 1) = Bz(0, Bz.get_Ny() - 1); // правый нижний узел


  //  for (i = 0ull; i < Nx; ++i)
  //    for (j = 0ull; j < Ny; ++j)
  //      for (k = 0ull; k < Nz; ++k)
  //      {
  //        Ex(i, j, k) = Ex(i, j, k) + C * E_dt * ((Bz(i, j, k) - Bz(i, j - 1, k)) / dy /*- (By(i, j, k) - By(i, j, k - 1)) / dz*/);
  //        Ey(i, j, k) = Ey(i, j, k) + C * E_dt * (/*(Bx(i, j, k) - Bx(i, j, k - 1)) / dz*/ - (Bz(i, j, k) - Bz(i - 1, j, k)) / dx);
  //        Ez(i, j, k) = Ez(i, j, k) + C * E_dt * ((By(i, j, k) - By(i - 1, j, k)) / dx - (Bx(i, j, k) - Bx(i, j - 1, k)) / dy);
  //      }
  //  // Обновляю границу - верх 2-х мерной системы
  //  for (int64_t x = 0; x < Ex.get_Nx(); ++x)
  //  {
  //    *(Ex.data() + x + 1) = Ex(x, Ex.get_Ny() - 1); // снизу
  //    Ex(x, Ex.get_Ny()) = Ex(x, 0); // сверху

  //    *(Ey.data() + x + 1) = Ey(x, Ey.get_Ny() - 1); // снизу
  //    Ey(x, Ey.get_Ny()) = Ey(x, 0); // сверху

  //    *(Ez.data() + x + 1) = Ez(x, Ez.get_Ny() - 1); // снизу
  //    Ez(x, Ez.get_Ny()) = Ez(x, 0); // сверху
  //  }
  //  for (int64_t y = 0; y < Ex.get_Ny(); ++y)
  //  {
  //    Ex(-1, y) = Ex(Ex.get_Nx() - 1, y); // слева
  //    Ex(Ex.get_Nx(), y) = Ex(0, y); // справа

  //    Ey(-1, y) = Ey(Ey.get_Nx() - 1, y); // слева
  //    Ey(Ey.get_Nx(), y) = Ey(0, y); // справа

  //    Ez(-1, y) = Ez(Ez.get_Nx() - 1, y); // слева
  //    Ez(Ez.get_Nx(), y) = Ez(0, y); // справа
  //  }
  //  Ex(-1, -1) = Ex(Ex.get_Nx() - 1, Ex.get_Ny() - 1); // левый нижний узел
  //  Ex(-1, Ex.get_Ny()) = Ex(Ex.get_Nx() - 1, 0); // левый верхний узел
  //  Ex(Ex.get_Nx(), Ex.get_Ny()) = Ex(0, 0); // правый верхний узел
  //  *(Ex.data() + Ex.get_Nx() + 1) = Ex(0, Ex.get_Ny() - 1); // правый нижний узел

  //  Ey(-1, -1) = Ey(Ey.get_Nx() - 1, Ey.get_Ny() - 1); // левый нижний узел
  //  Ey(-1, Ey.get_Ny()) = Ey(Ey.get_Nx() - 1, 0); // левый верхний узел
  //  Ey(Ey.get_Nx(), Ey.get_Ny()) = Ey(0, 0); // правый верхний узел
  //  *(Ey.data() + Ey.get_Nx() + 1) = Ey(0, Ey.get_Ny() - 1); // правый нижний узел

  //  Ez(-1, -1) = Ez(Ez.get_Nx() - 1, Ez.get_Ny() - 1); // левый нижний узел
  //  Ez(-1, Ez.get_Ny()) = Ez(Ez.get_Nx() - 1, 0); // левый верхний узел
  //  Ez(Ez.get_Nx(), Ez.get_Ny()) = Ez(0, 0); // правый верхний узел
  //  *(Ez.data() + Ez.get_Nx() + 1) = Ez(0, Ez.get_Ny() - 1); // правый нижний узел


  //  for (i = 0ull; i < Nx; ++i)
  //    for (j = 0ull; j < Ny; ++j)
  //      for (k = 0ull; k < Nz; ++k)
  //      {
  //        Bx(i, j, k) = Bx(i, j, k) + C * B_dt * (/*(Ey(i, j, k + 1) - Ey(i, j, k)) / dz*/ - (Ez(i, j + 1, k) - Ez(i, j, k)) / dy);
  //        By(i, j, k) = By(i, j, k) + C * B_dt * ((Ez(i + 1, j, k) - Ez(i, j, k)) / dx /*- (Ex(i, j, k + 1) - Ex(i, j, k)) / dz*/);
  //        Bz(i, j, k) = Bz(i, j, k) + C * B_dt * ((Ex(i, j + 1, k) - Ex(i, j, k)) / dy - (Ey(i + 1, j, k) - Ey(i, j, k)) / dx);
  //      }
  //  // Обновляю границу - верх 2-х мерной системы
  //  for (int64_t x = 0; x < Bx.get_Nx(); ++x)
  //  {
  //    *(Bx.data() + x + 1) = Bx(x, Bx.get_Ny() - 1); // снизу
  //    Bx(x, Bx.get_Ny()) = Bx(x, 0); // сверху

  //    *(By.data() + x + 1) = By(x, By.get_Ny() - 1); // снизу
  //    By(x, By.get_Ny()) = By(x, 0); // сверху

  //    *(Bz.data() + x + 1) = Bz(x, Bz.get_Ny() - 1); // снизу
  //    Bz(x, Bz.get_Ny()) = Bz(x, 0); // сверху
  //  }
  //  for (int64_t y = 0; y < Bx.get_Ny(); ++y)
  //  {
  //    Bx(-1, y) = Bx(Bx.get_Nx() - 1, y); // слева
  //    Bx(Bx.get_Nx(), y) = Bx(0, y); // справа

  //    By(-1, y) = By(By.get_Nx() - 1, y); // слева
  //    By(By.get_Nx(), y) = By(0, y); // справа

  //    Bz(-1, y) = Bz(Bz.get_Nx() - 1, y); // слева
  //    Bz(Bz.get_Nx(), y) = Bz(0, y); // справа
  //  }
  //  Bx(-1, -1) = Bx(Bx.get_Nx() - 1, Bx.get_Ny() - 1); // левый нижний узел
  //  Bx(-1, Bx.get_Ny()) = Bx(Bx.get_Nx() - 1, 0); // левый верхний узел
  //  Bx(Bx.get_Nx(), Bx.get_Ny()) = Bx(0, 0); // правый верхний узел
  //  *(Bx.data() + Bx.get_Nx() + 1) = Bx(0, Bx.get_Ny() - 1); // правый нижний узел

  //  By(-1, -1) = By(By.get_Nx() - 1, By.get_Ny() - 1); // левый нижний узел
  //  By(-1, By.get_Ny()) = By(By.get_Nx() - 1, 0); // левый верхний узел
  //  By(By.get_Nx(), By.get_Ny()) = By(0, 0); // правый верхний узел
  //  *(By.data() + By.get_Nx() + 1) = By(0, By.get_Ny() - 1); // правый нижний узел

  //  Bz(-1, -1) = Bz(Bz.get_Nx() - 1, Bz.get_Ny() - 1); // левый нижний узел
  //  Bz(-1, Bz.get_Ny()) = Bz(Bz.get_Nx() - 1, 0); // левый верхний узел
  //  Bz(Bz.get_Nx(), Bz.get_Ny()) = Bz(0, 0); // правый верхний узел
  //  *(Bz.data() + Bz.get_Nx() + 1) = Bz(0, Bz.get_Ny() - 1); // правый нижний узел
  //}

//#endif // MPI
}

void FDTD::FDTD::write_fields_to_file(const char* path, const Component E, const Component B, const double delta, const int64_t row_number)
{

  Axis axis = get_axis(E, B);
  switch (axis)
  {
  case Axis::Ox:
    Ex.write_field_to_file_OX(path, row_number);
    Ey.write_field_to_file_OX(path, row_number);
    Ez.write_field_to_file_OX(path, row_number);
    Bx.write_field_to_file_OX(path, row_number);
    By.write_field_to_file_OX(path, row_number);
    Bz.write_field_to_file_OX(path, row_number);
    break;
  case Axis::Oy:
    Ex.write_field_to_file_OY(path, row_number);
    Ey.write_field_to_file_OY(path, row_number);
    Ez.write_field_to_file_OY(path, row_number);
    Bx.write_field_to_file_OY(path, row_number);
    By.write_field_to_file_OY(path, row_number);
    Bz.write_field_to_file_OY(path, row_number);
    break;
  case Axis::Oz:
    Ex.write_field_to_file_OZ(path, row_number);
    Ey.write_field_to_file_OZ(path, row_number);
    Ez.write_field_to_file_OZ(path, row_number);
    Bx.write_field_to_file_OZ(path, row_number);
    By.write_field_to_file_OZ(path, row_number);
    Bz.write_field_to_file_OZ(path, row_number);
    break;
  default:
    break;
  }

  std::ofstream outfile;
  outfile.open(path, std::ios::app);
  if (!outfile.is_open())
  {
    std::cout << "The file can't be opened!" << std::endl;
    exit(-1);
  }
  outfile << delta << std::endl;
  outfile.close();
}

FDTD::Axis FDTD::FDTD::get_axis(const Component E, const Component B)
{
  if ((E == Component::Ey && B == Component::Bz) || (E == Component::Ez && B == Component::By)) return Axis::Ox;
  else if ((E == Component::Ez && B == Component::Bx) || (E == Component::Ex && B == Component::Bz)) return Axis::Oy;
  else if ((E == Component::Ex && B == Component::By) || (E == Component::Ey && B == Component::Bx)) return Axis::Oz;
  else
  {
    std::cout << "\nGet_Axis: Error! Wrong Components!\n";
    exit(-1);
  }
}
