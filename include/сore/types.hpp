#pragma once
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>


namespace FDTD {

constexpr double PI = 3.14159265358979323846;
constexpr double C = 29979245800.0; // speed of light (CGS)

namespace fs = std::filesystem;

constexpr const char *path_to_analytical_data_directory =
    PATH_TO_ANALYTICAL_DATA_DIRECTORY;
constexpr const char *path_to_calculated_data_directory =
    PATH_TO_CALCULATED_DATA_DIRECTORY;
constexpr const char *path_to_convergence_data = PATH_TO_CONVERGENCE_DATA;

// ---MPI---
constexpr int MPI_DIMENSION = 3;
// ---MPI---

using Field = std::vector<double>;
using TimeStep = double;
using UnboundedGridSizes = GridSizes;
using BoundedGridSizes = GridSizes;

enum class Axis { X, Y, Z };
enum class Component { Ex, Ey, Ez, Bx, By, Bz };
std::ostream &operator<<(std::ostream &os, const Component &comp);

// ------------------------------------------------------
enum class Shift { shifted, unshifted };
// ------------------------------------------------------

struct GridSizes {
  int Nx, Ny, Nz;
};

struct GridCoordinatesBounds {
  double ax, ay, az;
  double bx, by, bz;
};

struct GridCoordinatesSteps {
  double dx, dy, dz;
};

} // namespace FDTD
