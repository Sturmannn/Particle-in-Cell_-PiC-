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
  Bz(_fields.By) {
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

uint64_t FDTD::FDTD::get_Nx(void) const { return Nx; }

uint64_t FDTD::FDTD::get_Ny(void) const { return Ny; }

void FDTD::FDTD::field_update(const double t) {
  for (double time = 0.0; time < t; time += dt)
    for (uint64_t j = 0ll; j < Ny; ++j)
      for (uint64_t i = 0ll; i < Nx; ++i) {
        Ex(i, j) =
          Ex(i, j) + C * dt * 0.5 * (Bz(i, j + 1ll) - Bz(i, j - 1ll)) / dy;
        Ey(i, j) =
          Ey(i, j) - C * dt * 0.5 * (Bz(i + 1ll, j) - Bz(i - 1ll, j)) / dx;
        Ez(i, j) = Ez(i, j) + C * dt * 0.5 *
          ((By(i + 1ll, j) - By(i - 1ll, j)) / dx -
            (Bx(i, j + 1ll) - Bx(i, j - 1ll)) / dy);

        Bx(i, j) =
          Bx(i, j) - C * dt * 0.5 * (Ez(i, j + 1ll) - Ez(i, j - 1ll)) / dy;
        By(i, j) =
          By(i, j) + C * dt * 0.5 * (Ez(i + 1ll, j) - Ez(i - 1ll, j)) / dx;
        Bz(i, j) = Bz(i, j) - C * dt * 0.5 *
          ((Ey(i + 1ll, j) - Ey(i - 1ll, j)) / dx -
            (Ex(i, j + 1ll) - Ex(i, j - 1ll)) / dy);
      }
}
