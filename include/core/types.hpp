#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstring>


namespace FDTD {

constexpr double PI = 3.14159265358979323846;
constexpr double C = 29979245800.0; // speed of light (CGS)

namespace fs = std::filesystem;

constexpr const char *path_to_analytical_data_directory =
    PATH_TO_ANALYTICAL_DATA_DIRECTORY;
constexpr const char *path_to_calculated_data_directory =
    PATH_TO_CALCULATED_DATA_DIRECTORY;
constexpr const char *path_to_convergence_data = PATH_TO_CONVERGENCE_DATA;
constexpr const char *path_to_measurements_file = PATH_TO_MEASUREMENTS_FILE;

// ---MPI---
constexpr int MPI_DIMENSION = 3;
// ---MPI---

enum class Axis { X, Y, Z };
enum class Component { Ex, Ey, Ez, Bx, By, Bz };
inline std::ostream &operator<<(std::ostream &os, const Component &comp) {
  switch(comp) {
    case Component::Ex: os << "Ex"; break;
    case Component::Ey: os << "Ey"; break;
    case Component::Ez: os << "Ez"; break;
    case Component::Bx: os << "Bx"; break;
    case Component::By: os << "By"; break;
    case Component::Bz: os << "Bz"; break;
    default: os << "Unknown"; break;
}
return os;
};

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

using Field = std::vector<double>;
using TimeStep = double;
using UnboundedGridSizes = GridSizes;
using BoundedGridSizes = GridSizes;

} // namespace FDTD
