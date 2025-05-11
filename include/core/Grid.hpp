#pragma once

#include "types.hpp"

namespace FDTD {

class Grid {
public:
  explicit Grid(const UnboundedGridSizes &_unbounded_sizes,
                const GridCoordinatesBounds &_bounds, const TimeStep &_dt)
      : unbounded_grid_sizes(_unbounded_sizes), bounds(_bounds), dt(_dt) {

    coordinatesSteps.dx =
        (bounds.bx - bounds.ax) / static_cast<double>(unbounded_grid_sizes.Nx);
    coordinatesSteps.dy =
        (bounds.by - bounds.ay) / static_cast<double>(unbounded_grid_sizes.Ny);
    coordinatesSteps.dz =
        (bounds.bz - bounds.az) / static_cast<double>(unbounded_grid_sizes.Nz);

    bounded_grid_sizes.Nx = unbounded_grid_sizes.Nx + 2;
    bounded_grid_sizes.Ny = unbounded_grid_sizes.Ny + 2;
    bounded_grid_sizes.Nz = unbounded_grid_sizes.Nz + 2;
  };

  Grid(const Grid &other) = default;

  UnboundedGridSizes get_unbounded_grid_sizes() const {
    return unbounded_grid_sizes;
  }
  UnboundedGridSizes get_mpi_local_unbounded_subdomain_sizes() const {
    return mpi_local_unbounded_subdomain_sizes;
  }
  int mpi_get_total_bounded_subdomain_sizes() const {
    return mpi_local_bounded_subdomain_sizes.Nx *
           mpi_local_bounded_subdomain_sizes.Ny *
           mpi_local_bounded_subdomain_sizes.Nz;
  }
  void set_mpi_local_bounded_subdomain_sizes(
      const BoundedGridSizes &mpi_local_bounded_subdomain_sizes) {
        this->mpi_local_bounded_subdomain_sizes.Nx = mpi_local_bounded_subdomain_sizes.Nx;
        this->mpi_local_bounded_subdomain_sizes.Ny = mpi_local_bounded_subdomain_sizes.Ny;
        this->mpi_local_bounded_subdomain_sizes.Nz = mpi_local_bounded_subdomain_sizes.Nz;
  }
  void set_mpi_local_unbounded_subdomain_sizes(
      const UnboundedGridSizes &mpi_local_unbounded_subdomain_sizes) {
    this->mpi_local_unbounded_subdomain_sizes.Nx = mpi_local_unbounded_subdomain_sizes.Nx;
    this->mpi_local_unbounded_subdomain_sizes.Ny = mpi_local_unbounded_subdomain_sizes.Ny;
    this->mpi_local_unbounded_subdomain_sizes.Nz = mpi_local_unbounded_subdomain_sizes.Nz;
  }

  BoundedGridSizes get_bounded_grid_sizes() const { return bounded_grid_sizes; }

  BoundedGridSizes get_mpi_local_bounded_subdomain_sizes() const {
    return mpi_local_bounded_subdomain_sizes;
  }

  GridCoordinatesBounds get_bounds() const { return bounds; }
  GridCoordinatesSteps get_coordinates_steps() const {
    return coordinatesSteps;
  }
  TimeStep get_dt() const { return dt; }

  double get_index(const int i, const int j, const int k) {
    int Nx = mpi_local_bounded_subdomain_sizes.Nx;
    int Ny = mpi_local_bounded_subdomain_sizes.Ny;
    int Nz = mpi_local_bounded_subdomain_sizes.Nz;

    return Ny * Nz * (i + 1) + Nz * (j + 1) + (k + 1);
  }

private:
  BoundedGridSizes bounded_grid_sizes;     // Nx + 2, Ny + 2, Nz + 2
  UnboundedGridSizes unbounded_grid_sizes; // Nx, Ny, Nz
  // ---MPI---
  BoundedGridSizes mpi_local_bounded_subdomain_sizes; // Nx + 2, Ny + 2, Nz + 2
  UnboundedGridSizes mpi_local_unbounded_subdomain_sizes; // Nx, Ny, Nz
  // ---MPI---
  GridCoordinatesBounds bounds;          // ax, ay, az, bx, by, bz
  GridCoordinatesSteps coordinatesSteps; // dx, dy, dz
  TimeStep dt;                           // dt
};

} // namespace FDTD