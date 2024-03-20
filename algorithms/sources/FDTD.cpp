#include "FDTD.hpp"

FDTD::FDTD::FDTD(const std::tuple<uint64_t, uint64_t, uint64_t>& Nx_Ny_Nz,
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

void FDTD::FDTD::set_Nx_Ny_Nz(uint64_t _Nx, uint64_t _Ny, uint64_t _Nz)
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

  uint64_t i = 0ull;
  uint64_t j = 0ull;

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

void FDTD::FDTD::field_update(const uint64_t t)
{
  if (dt == 0.0)
  {
    std::cout << "Time step is null";
    exit(-1);
  }

  uint64_t i = 0ull;
  uint64_t j = 0ull;

  for (uint64_t time = 0ull; time < t; time++)
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
  uint64_t i = 0ull;
  uint64_t j = 0ull;
  uint64_t k = 0ull;

  //double start = omp_get_wtime();
  for (double time = 0.0; time < t; time += dt) // ��������� �� time += E_dt
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

void FDTD::FDTD::shifted_field_update(const uint64_t t)
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
#ifdef MPI

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  uint64_t Nx_local = get_Nx();
  uint64_t Ny_local = get_Ny() / size;

  Field::ComputingField Ex_local(Nx_local, Ny_local);
  Field::ComputingField Ey_local(Nx_local, Ny_local);
  Field::ComputingField Ez_local(Nx_local, Ny_local);
  Field::ComputingField Bx_local(Nx_local, Ny_local);
  Field::ComputingField By_local(Nx_local, Ny_local);
  Field::ComputingField Bz_local(Nx_local, Ny_local);

  // ���� ��� ��������������, ��� Ny ������ ���������� ��������� size! ����� ���� ��� ��������������� MPI ��� ���������� ������!
  MPI_Scatter(Ex.data(), Nx_local * Ny_local, MPI_DOUBLE, Ex_local.data(), Nx_local * Ny_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(Ey.data(), Nx_local * Ny_local, MPI_DOUBLE, Ey_local.data(), Nx_local * Ny_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(Ez.data(), Nx_local * Ny_local, MPI_DOUBLE, Ez_local.data(), Nx_local * Ny_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Scatter(Bx.data(), Nx_local * Ny_local, MPI_DOUBLE, Bx_local.data(), Nx_local * Ny_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(By.data(), Nx_local * Ny_local, MPI_DOUBLE, By_local.data(), Nx_local * Ny_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(Bz.data(), Nx_local * Ny_local, MPI_DOUBLE, Bz_local.data(), Nx_local * Ny_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for (uint64_t time = 0ull; time < t; time++)
  {
    for (i = 0; i < Nx; ++i)
      for (j = 0; j < Ny; ++j)
        for (k = 0; k < Nz; ++k)
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


#else  
  for (uint64_t time = 0ull; time < t; time++)
  {
    for (i = 0; i < Nx; ++i)
      for (j = 0; j < Ny; ++j)
        for (k = 0; k < Nz; ++k)
        {
          Bx(i, j, k) = Bx(i, j, k) + C * B_dt * (/*(Ey(i, j, k + 1) - Ey(i, j, k)) / dz*/ - (Ez(i, j + 1, k) - Ez(i, j, k)) / dy);
          By(i, j, k) = By(i, j, k) + C * B_dt * ((Ez(i + 1, j, k) - Ez(i, j, k)) / dx /*- (Ex(i, j, k + 1) - Ex(i, j, k)) / dz*/);
          Bz(i, j, k) = Bz(i, j, k) + C * B_dt * ((Ex(i, j + 1, k) - Ex(i, j, k)) / dy - (Ey(i + 1, j, k) - Ey(i, j, k)) / dx);
          //std::cout << "i = " << i << " j = " << j << " k = " << k << '\n';
        }
    // �������� ������� - ���� 2-� ������ �������
    for (uint64_t x = 0; x < Bx.get_Nx(); ++x)
    {
      *(Bx.data() + x + 1) = Bx(x, Bx.get_Ny() - 1); // �����
      Bx(x, Bx.get_Ny()) = Bx(x, 0); // ������

      *(By.data() + x + 1) = By(x, By.get_Ny() - 1); // �����
      By(x, By.get_Ny() - 1) = By(x, 0); // ������

      *(Bz.data() + x + 1) = Bz(x, Bz.get_Ny() - 1); // �����
      Bz(x, Bz.get_Ny()) = Bz(x, 0); // ������
    }
    for (uint64_t y = 0; y < Bx.get_Ny(); ++y)
    {
      Bx(-1, y) = Bx(Bx.get_Nx() - 1, y); // �����
      Bx(Bx.get_Nx(), y) = Bx(0, y); // ������

      By(-1, y) = By(By.get_Nx() - 1, y); // �����
      By(By.get_Nx(), y) = By(0, y); // ������

      Bz(-1, y) = Bz(Bz.get_Nx() - 1, y); // �����
      Bz(Bz.get_Nx(), y) = Bz(0, y); // ������
    }
    Bx(-1, -1) = Bx(Bx.get_Nx() - 1, Bx.get_Ny() - 1); // ����� ������ ����
    Bx(-1, Bx.get_Ny()) = Bx(Bx.get_Nx() - 1, 0); // ����� ������� ����
    Bx(Bx.get_Nx(), Bx.get_Ny()) = Bx(0, 0); // ������ ������� ����
    *(Bx.data() + Bx.get_Nx() + 1) = Bx(0, Bx.get_Ny() - 1); // ������ ������ ����

    By(-1, -1) = By(By.get_Nx() - 1, By.get_Ny() - 1); // ����� ������ ����
    By(-1, By.get_Ny()) = By(By.get_Nx() - 1, 0); // ����� ������� ����
    By(By.get_Nx(), By.get_Ny()) = By(0, 0); // ������ ������� ����
    *(By.data() + By.get_Nx() + 1) = By(0, By.get_Ny() - 1); // ������ ������ ����

    Bz(-1, -1) = Bz(Bz.get_Nx() - 1, Bz.get_Ny() - 1); // ����� ������ ����
    Bz(-1, Bz.get_Ny()) = Bz(Bz.get_Nx() - 1, 0); // ����� ������� ����
    Bz(Bz.get_Nx(), Bz.get_Ny()) = Bz(0, 0); // ������ ������� ����
    *(Bz.data() + Bz.get_Nx() + 1) = Bz(0, Bz.get_Ny() - 1); // ������ ������ ����


    for (i = 0ull; i < Nx; ++i)
      for (j = 0ull; j < Ny; ++j)
        for (k = 0ull; k < Nz; ++k)
        {
          Ex(i, j, k) = Ex(i, j, k) + C * E_dt * ((Bz(i, j, k) - Bz(i, j - 1, k)) / dy /*- (By(i, j, k) - By(i, j, k - 1)) / dz*/);
          Ey(i, j, k) = Ey(i, j, k) + C * E_dt * (/*(Bx(i, j, k) - Bx(i, j, k - 1)) / dz*/ - (Bz(i, j, k) - Bz(i - 1, j, k)) / dx);
          Ez(i, j, k) = Ez(i, j, k) + C * E_dt * ((By(i, j, k) - By(i - 1, j, k)) / dx - (Bx(i, j, k) - Bx(i, j - 1, k)) / dy);
        }
    // �������� ������� - ���� 2-� ������ �������
    for (uint64_t x = 0; x < Ex.get_Nx(); ++x)
    {
      *(Ex.data() + x + 1) = Ex(x, Ex.get_Ny() - 1); // �����
      Ex(x, Ex.get_Ny()) = Ex(x, 0); // ������

      *(Ey.data() + x + 1) = Ey(x, Ey.get_Ny() - 1); // �����
      Ey(x, Ey.get_Ny()) = Ey(x, 0); // ������

      *(Ez.data() + x + 1) = Ez(x, Ez.get_Ny() - 1); // �����
      Ez(x, Ez.get_Ny()) = Ez(x, 0); // ������
    }
    for (uint64_t y = 0; y < Ex.get_Ny(); ++y)
    {
      Ex(-1, y) = Ex(Ex.get_Nx() - 1, y); // �����
      Ex(Ex.get_Nx(), y) = Ex(0, y); // ������

      Ey(-1, y) = Ey(Ey.get_Nx() - 1, y); // �����
      Ey(Ey.get_Nx(), y) = Ey(0, y); // ������

      Ez(-1, y) = Ez(Ez.get_Nx() - 1, y); // �����
      Ez(Ez.get_Nx(), y) = Ez(0, y); // ������
    }
    Ex(-1, -1) = Ex(Ex.get_Nx() - 1, Ex.get_Ny() - 1); // ����� ������ ����
    Ex(-1, Ex.get_Ny()) = Ex(Ex.get_Nx() - 1, 0); // ����� ������� ����
    Ex(Ex.get_Nx(), Ex.get_Ny()) = Ex(0, 0); // ������ ������� ����
    *(Ex.data() + Ex.get_Nx() + 1) = Ex(0, Ex.get_Ny() - 1); // ������ ������ ����

    Ey(-1, -1) = Ey(Ey.get_Nx() - 1, Ey.get_Ny() - 1); // ����� ������ ����
    Ey(-1, Ey.get_Ny()) = Ey(Ey.get_Nx() - 1, 0); // ����� ������� ����
    Ey(Ey.get_Nx(), Ey.get_Ny()) = Ey(0, 0); // ������ ������� ����
    *(Ey.data() + Ey.get_Nx() + 1) = Ey(0, Ey.get_Ny() - 1); // ������ ������ ����

    Ez(-1, -1) = Ez(Ez.get_Nx() - 1, Ez.get_Ny() - 1); // ����� ������ ����
    Ez(-1, Ez.get_Ny()) = Ez(Ez.get_Nx() - 1, 0); // ����� ������� ����
    Ez(Ez.get_Nx(), Ez.get_Ny()) = Ez(0, 0); // ������ ������� ����
    *(Ez.data() + Ez.get_Nx() + 1) = Ez(0, Ez.get_Ny() - 1); // ������ ������ ����


    for (i = 0ull; i < Nx; ++i)
      for (j = 0ull; j < Ny; ++j)
        for (k = 0ull; k < Nz; ++k)
        {
          Bx(i, j, k) = Bx(i, j, k) + C * B_dt * (/*(Ey(i, j, k + 1) - Ey(i, j, k)) / dz*/ - (Ez(i, j + 1, k) - Ez(i, j, k)) / dy);
          By(i, j, k) = By(i, j, k) + C * B_dt * ((Ez(i + 1, j, k) - Ez(i, j, k)) / dx /*- (Ex(i, j, k + 1) - Ex(i, j, k)) / dz*/);
          Bz(i, j, k) = Bz(i, j, k) + C * B_dt * ((Ex(i, j + 1, k) - Ex(i, j, k)) / dy - (Ey(i + 1, j, k) - Ey(i, j, k)) / dx);
        }
    // �������� ������� - ���� 2-� ������ �������
    for (uint64_t x = 0; x < Bx.get_Nx(); ++x)
    {
      *(Bx.data() + x + 1) = Bx(x, Bx.get_Ny() - 1); // �����
      Bx(x, Bx.get_Ny()) = Bx(x, 0); // ������

      *(By.data() + x + 1) = By(x, By.get_Ny() - 1); // �����
      By(x, By.get_Ny()) = By(x, 0); // ������

      *(Bz.data() + x + 1) = Bz(x, Bz.get_Ny() - 1); // �����
      Bz(x, Bz.get_Ny()) = Bz(x, 0); // ������
    }
    for (uint64_t y = 0; y < Bx.get_Ny(); ++y)
    {
      Bx(-1, y) = Bx(Bx.get_Nx() - 1, y); // �����
      Bx(Bx.get_Nx(), y) = Bx(0, y); // ������

      By(-1, y) = By(By.get_Nx() - 1, y); // �����
      By(By.get_Nx(), y) = By(0, y); // ������

      Bz(-1, y) = Bz(Bz.get_Nx() - 1, y); // �����
      Bz(Bz.get_Nx(), y) = Bz(0, y); // ������
    }
    Bx(-1, -1) = Bx(Bx.get_Nx() - 1, Bx.get_Ny() - 1); // ����� ������ ����
    Bx(-1, Bx.get_Ny()) = Bx(Bx.get_Nx() - 1, 0); // ����� ������� ����
    Bx(Bx.get_Nx(), Bx.get_Ny()) = Bx(0, 0); // ������ ������� ����
    *(Bx.data() + Bx.get_Nx() + 1) = Bx(0, Bx.get_Ny() - 1); // ������ ������ ����

    By(-1, -1) = By(By.get_Nx() - 1, By.get_Ny() - 1); // ����� ������ ����
    By(-1, By.get_Ny()) = By(By.get_Nx() - 1, 0); // ����� ������� ����
    By(By.get_Nx(), By.get_Ny()) = By(0, 0); // ������ ������� ����
    *(By.data() + By.get_Nx() + 1) = By(0, By.get_Ny() - 1); // ������ ������ ����

    Bz(-1, -1) = Bz(Bz.get_Nx() - 1, Bz.get_Ny() - 1); // ����� ������ ����
    Bz(-1, Bz.get_Ny()) = Bz(Bz.get_Nx() - 1, 0); // ����� ������� ����
    Bz(Bz.get_Nx(), Bz.get_Ny()) = Bz(0, 0); // ������ ������� ����
    *(Bz.data() + Bz.get_Nx() + 1) = Bz(0, Bz.get_Ny() - 1); // ������ ������ ����
  }
#endif // MPI
}

void FDTD::FDTD::write_fields_to_file(const char* path, const Component E, const Component B, const double delta, const uint64_t row_number)
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
