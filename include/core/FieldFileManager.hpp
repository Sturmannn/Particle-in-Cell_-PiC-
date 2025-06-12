#pragma once

#include "FDTD_MPI.hpp"
#include "types.hpp"

namespace FDTD {

class FieldFileManager {
public:
  FieldFileManager(std::shared_ptr<MPI_Wrapper> _mpi_wrapper) : mpi_wrapper(_mpi_wrapper) {}
  void clear_files(const char *directory_path);

  void write_fields_to_file(const char *directory_path, Component E,
                            Component B, const double delta,
                            const FDTD_MPI &fields,
                            const int row_number = 0); // The col is fixed

private:
    std::shared_ptr<MPI_Wrapper> mpi_wrapper;
};

inline void FieldFileManager::clear_files(const char *directory_path) {

  int rank = mpi_wrapper->get_world_rank();
  int size = mpi_wrapper->get_world_size();

  fs::path dir_path = fs::path(directory_path);
  fs::path file_path = fs::path(directory_path);
  
  if (file_path.parent_path().filename() == "my_data")
  file_path /= "my_data_" + std::to_string(rank) + ".csv";
  else if (file_path.parent_path().filename() == "analytical_data")
  file_path /= "analytical_data_" + std::to_string(rank) + ".csv";
  else {
    std::cout
    << "Field::ComputingField::clear_files() error: invalid directory path"
    << std::endl;
    exit(-1);
  }
  
  // Clearing or Creating the file
  try {
    std::ofstream file(file_path, std::ofstream::out | std::ofstream::trunc);
    if (file.is_open()) {
      file.close();
    } else {
      std::cout << "Process " << rank << " failed to open file: " << file_path
                << std::endl;
    }
  } catch (const std::exception &e) {
    std::cout << "Error in process " << rank << " while handling file "
              << file_path << ": " << e.what() << std::endl;
  }

  // Removing extra files
  if (rank == 0) {
    int i = size;
    while (true) {
      fs::path extra_file_path;
      if (dir_path.parent_path().filename() == "my_data")
        extra_file_path = dir_path / ("my_data_" + std::to_string(i) + ".csv");
      else if (dir_path.parent_path().filename() == "analytical_data")
        extra_file_path =
            dir_path / ("analytical_data_" + std::to_string(i) + ".csv");

      if (fs::exists(extra_file_path)) {
        try {
          fs::remove(extra_file_path);
        } catch (const std::exception &e) {
          std::cerr << "Error while removing extra file " << extra_file_path
                    << ": " << e.what() << std::endl;
        }
      } else {
        break;
      }
      ++i;
    }
  }
}

inline void FieldFileManager::write_fields_to_file(const char *directory_path,
                                                   Component E, Component B,
                                                   const double delta,
                                                   const FDTD_MPI &fields,
                                                   const int row_number) {
  int rank = mpi_wrapper->get_world_rank();

  fs::path file_path = fs::path(directory_path);

  if (file_path.parent_path().filename() == "my_data")
    file_path /= "my_data_" + std::to_string(rank) + ".csv";
  else if (file_path.parent_path().filename() == "analytical_data")
    file_path /= "analytical_data_" + std::to_string(rank) + ".csv";
  else {
    std::cout
        << "FDTD::FDTD::write_fields_to_file() error: invalid directory path"
        << std::endl;
    exit(-1);
  }

  std::ofstream outfile;
  outfile.open(file_path, std::ios::app);

  Axis axis = fields.get_axis(E, B);
  int index = 0;
  UnboundedGridSizes grid_sizes = fields.get_grid()->get_mpi_local_unbounded_subdomain_sizes();
  int Nx = grid_sizes.Nx;
  int Ny = grid_sizes.Ny;
  int Nz = grid_sizes.Nz;
  std::vector<Field> fields_vector = {fields.EX(), fields.EY(), fields.EZ(),
                                      fields.BX(), fields.BY(), fields.BZ()};
  switch (axis) {
  case Axis::X:
    if (row_number >= Nx) {
      std::cout << "The row number is out of range!" << std::endl;
      exit(-1);
    }
    for (const auto field : fields_vector) {
      for (int i = 0; i < Nx; ++i) {
        index = fields.get_grid()->get_index(i, row_number, 0);
        outfile << field[index] << (i < Nx - 1 ? ';' : '\n');
      }
    }
    outfile << delta << std::endl;
    break;
  case Axis::Y:
    if (row_number >= Ny) {
      std::cout << "The row number is out of range!" << std::endl;
      exit(-1);
    }
    for (const auto field : fields_vector) {
      for (int j = 0; j < Ny; ++j) {
        index = fields.get_grid()->get_index(row_number, j, 0);
        outfile << field[index] << (j < Ny - 1 ? ';' : '\n');
      }
    }
    outfile << delta << std::endl;
    break;
  case Axis::Z:
    if (row_number >= Nz) {
      std::cout << "The row number is out of range!" << std::endl;
      exit(-1);
    }
    for (const auto field : fields_vector) {
      for (int k = 0; k < Nz; ++k) {
        index = fields.get_grid()->get_index(0, row_number, k);
        outfile << field[index] << (k < Nz - 1 ? ';' : '\n');
      }
    }
    outfile << delta << std::endl;
    break;
  default:
    break;
  }
  outfile.close();
}

} // namespace FDTD