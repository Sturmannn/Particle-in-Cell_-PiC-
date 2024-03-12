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
  for (double time = 0.0; time < t; time += dt) // ÏÎÏÐÀÂÈÒÜ ÍÀ time += E_dt
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
  uint64_t i = 0ull;
  uint64_t j = 0ull;
  uint64_t k = 0ull;
  //double start = omp_get_wtime();
  for (uint64_t time = 0ull; time < t; time++)
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
}

void FDTD::FDTD::write_fields_to_file_OX(const char* path, const double dx, uint64_t j)
{
  Ex.write_field_to_file_OX(path);
  Ey.write_field_to_file_OX(path);
  Ez.write_field_to_file_OX(path);
  Bx.write_field_to_file_OX(path);
  By.write_field_to_file_OX(path);
  Bz.write_field_to_file_OX(path);

  std::ofstream outfile;
  outfile.open(path, std::ios::app);
  if (!outfile.is_open())
  {
    std::cout << "The file can't be opened!" << std::endl;
    exit(-1);
  }
  outfile << dx << std::endl;
  outfile.close();
}

void FDTD::FDTD::write_fields_to_file_OY(const char* path, const double dy, uint64_t i)
{
  Ex.write_field_to_file_OY(path);
  Ey.write_field_to_file_OY(path);
  Ez.write_field_to_file_OY(path);
  Bx.write_field_to_file_OY(path);
  By.write_field_to_file_OY(path);
  Bz.write_field_to_file_OY(path);

  std::ofstream outfile;
  outfile.open(path, std::ios::app);
  if (!outfile.is_open())
  {
    std::cout << "The file can't be opened!" << std::endl;
    exit(-1);
  }
  outfile << dy << std::endl;
  outfile.close();
}

void FDTD::FDTD::write_fields_to_file_OZ(const char* path, const double dz, uint64_t k)
{
  Ex.write_field_to_file_OZ(path);
  Ey.write_field_to_file_OZ(path);
  Ez.write_field_to_file_OZ(path);
  Bx.write_field_to_file_OZ(path);
  By.write_field_to_file_OZ(path);
  Bz.write_field_to_file_OZ(path);

  std::ofstream outfile;
  outfile.open(path, std::ios::app);
  if (!outfile.is_open())
  {
    std::cout << "The file can't be opened!" << std::endl;
    exit(-1);
  }
  outfile << dz << std::endl;
  outfile.close();
}
