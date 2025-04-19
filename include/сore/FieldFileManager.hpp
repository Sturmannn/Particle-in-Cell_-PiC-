#pragma once

#include "FDTD_MPI.hpp"
#include "types.hpp"

namespace FDTD {

class FieldFileManager {
public:
  static void clear_files(const char *directory_path);

  void write_fields_to_file(const char *directory_path, Component E,
                            Component B, const double delta, const FDTD_MPI &fields,
                            const int64_t row_number = 0); // The col is fixed

private:
};

inline void FieldFileManager::clear_files(const char *directory_path) {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  fs::path dir_path = fs::path(directory_path);

  // Формирование имени файла для текущего процесса
  // Обращение сначала к parent_path из-за того, что путь заканчивается на /
  // (поэтому последний элемент пустой)
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

  // Очистка файла или его создание
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

  // Удаление лишних файлов
  // В случае ручного вмешательства в директории, следует очистить все файлы
  // вручную
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
        // Прекращение цикла, если файл не существует
        // Так как файлы идут по порядку, то если один не существует, то и
        // последующие тоже не существуют
        break;
      }
      ++i;
    }
  }
}

inline void FieldFileManager::write_fields_to_file(const char *directory_path,
                                                   Component E, Component B,
                                                   const double delta, const FDTD_MPI &fields,
                                                   const int64_t row_number) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  fs::path file_path = fs::path(directory_path);

  // Формирование имени файла для текущего процесса
  // Обращение сначала к patent_path из-за того, что путь заканчивается на /
  // (поэтому последний элемент пустой)
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
  UnboundedGridSizes grid_sizes = fields.get_grid()->get_unbounded_grid_sizes();
  int Nx = grid_sizes.Nx;
  int Ny = grid_sizes.Ny;
  int Nz = grid_sizes.Nz;
  std::vector<Field> fields_vector = {fields.EX(), fields.EY(), fields.EZ(),
                                       fields.BX(), fields.BY(), fields.BZ()};
  switch (axis) {
  case Axis::X:
    for (const auto field : fields_vector) {
        for (int i = 0; i < Nx; ++i){
            index = fields.get_grid()->get_index(i, index, 1);
            outfile << field[index] << (i < Nx - 1 ? ';' : '\n');
        }
    }
    break;
  case Axis::Y:
    for(const auto field : fields_vector) {
        for (int j = 0; j < Ny; ++j){
            index = fields.get_grid()->get_index(index, j, -1);
            outfile << field[index] << (j < Ny - 1 ? ';' : '\n');
        }
    }
    break;
  case Axis::Z:
    for (const auto field : fields_vector) {
        for (int k = 0; k < Nz; ++k){
            index = fields.get_grid()->get_index(0, index, k);
            outfile << field[index] << (k < Nz - 1 ? ';' : '\n');
        }
    }
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

} // namespace FDTD