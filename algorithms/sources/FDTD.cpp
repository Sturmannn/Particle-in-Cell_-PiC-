#include "FDTD.hpp"

FDTD::FDTD::FDTD(const std::pair<uint64_t, uint64_t>& Nx_Ny,
  const std::pair<double, double>& ax_ay,
  const std::pair<double, double>& bx_by, double _dt)
  : Ex(Nx_Ny.first, Nx_Ny.second),
  Ey(Nx_Ny.first, Nx_Ny.second),
  Ez(Nx_Ny.first, Nx_Ny.second),
  Bx(Nx_Ny.first, Nx_Ny.second),
  By(Nx_Ny.first, Nx_Ny.second),
  Bz(Nx_Ny.first, Nx_Ny.second) {
  Nx = Nx_Ny.first;
  Ny = Nx_Ny.second;

  ax = ax_ay.first;
  ay = ax_ay.second;

  bx = bx_by.first;
  by = bx_by.second;

  dx = (bx - ax) / Nx;
  dy = (by - ay) / Ny;

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

  ax = _fields.ax;
  ay = _fields.ay;

  bx = _fields.bx;
  by = _fields.by;

  dx = _fields.dx;
  dy = _fields.dy;

  dt = _fields.dt;
}

FDTD::FDTD::FDTD(FDTD&& _fields) noexcept
  : Ex(std::move(_fields.Ex)),
  Ey(std::move(_fields.Ey)),
  Ez(std::move(_fields.Ez)),
  Bx(std::move(_fields.Bx)),
  By(std::move(_fields.By)),
  Bz(std::move(_fields.Bz)) {
  
  Nx = std::move(_fields.Nx);
  Ny = std::move(_fields.Ny);

  ax = std::move(_fields.ax);
  ay = std::move(_fields.ay);
  
  bx = std::move(_fields.bx);
  by = std::move(_fields.by);
  
  dx = std::move(_fields.dx);
  dy = std::move(_fields.dy);
  dt = std::move(_fields.dt);
}

FDTD::FDTD& FDTD::FDTD::operator=(const FDTD& _fields)
{
  if (this != &_fields)
  {
    Ex = _fields.Ex;
    Ey = _fields.Ey;
    Ez = _fields.Ez;
    
    Bx = _fields.Bx;
    By = _fields.By;
    Bz = _fields.Bz;

    Nx = _fields.Nx;
    Ny = _fields.Ny;

    ax = _fields.ax;
    ay = _fields.ay;

    bx = _fields.bx;
    by = _fields.by;

    dx = _fields.dx;
    dy = _fields.dy;

    dt = _fields.dt;
  }
  return *this;
}

FDTD::FDTD& FDTD::FDTD::operator=(FDTD&& _fields) noexcept
{
  if (this != &_fields)
  {
    Ex = std::move(_fields.Ex);
    Ey = std::move(_fields.Ey);
    Ez = std::move(_fields.Ez);
    Bx = std::move(_fields.Bx);
    By = std::move(_fields.By);
    Bz = std::move(_fields.Bz);

    Nx = std::move(Nx);
    Ny = std::move(Ny);

    ax = std::move(ax);
    ay = std::move(ay);
    bx = std::move(bx);
    by = std::move(by);

    dx = std::move(dx);
    dy = std::move(dy);
    dt = std::move(dt);
  }
  return *this;
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
  //double start = omp_get_wtime();
  for (double time = 0.0; time < t; time += dt) // ÏÎÏÐÀÂÈÒÜ ÍÀ time += E_dt
  {
#pragma omp parallel for collapse(2)
    for (j = 0ull; j < Ny; ++j)
      for (i = 0ull; i < Nx; ++i)
      {
        Bx(i, j) = Bx(i, j) + C * B_dt * ((Ez(i, j) - Ez(i, j + 1ull)) / dy);
        By(i, j) = By(i, j) + C * B_dt * ((Ez(i + 1ull, j) - Ez(i, j)) / dx);
        Bz(i, j) = Bz(i, j) + C * B_dt * (((Ex(i, j + 1ull) - Ex(i, j)) / dy) - (Ey(i + 1ull, j) - Ey(i, j)) / dx);
      }
#pragma omp parallel for collapse(2)
    for (j = 0ull; j < Ny; ++j)
      for (i = 0ull; i < Nx; ++i)
      {
        Ex(i, j) = Ex(i, j) + C * E_dt * ((Bz(i, j) - Bz(i, j - 1ull)) / dy);
        Ey(i, j) = Ey(i, j) + C * E_dt * ((Bz(i - 1ull, j) - Bz(i, j)) / dx);
        Ez(i, j) = Ez(i, j) + C * E_dt * (((By(i, j) - By(i - 1ull, j)) / dx) - (Bx(i, j) - Bx(i, j - 1ull)) / dy);
      }
#pragma omp parallel for collapse(2)
    for (j = 0ull; j < Ny; ++j)
      for (i = 0ull; i < Nx; ++i)
      {
        Bx(i, j) = Bx(i, j) + C * B_dt * (Ez(i, j) - Ez(i, j + 1ull)) / dy;
        By(i, j) = By(i, j) + C * B_dt * (Ez(i + 1ull, j) - Ez(i, j)) / dx;
        Bz(i, j) = Bz(i, j) + C * B_dt * (((Ex(i, j + 1ull) - Ex(i, j)) / dy) - (Ey(i + 1ull, j) - Ey(i, j)) / dx);
      }
  }
  //double end = omp_get_wtime();
  //std::cout << "Time = " << end - start << std::endl;
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
