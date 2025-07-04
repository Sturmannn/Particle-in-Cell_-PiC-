#include "Field.hpp"

Field::ComputingField::ComputingField(const int64_t _Nx, const int64_t _Ny, const int64_t _Nz)
{
  if (_Nx <= 0 || _Ny <= 0 || _Nz < 0)
  {
    //throw std::runtime_error("Error: The grid size is incorrect");
    std::cout << "\nError: The grid size is incorrect\n";
    exit(-1);
  }
  Nx = _Nx;
  Ny = _Ny;
  _Nz == 0 ? Nz = 1 : Nz = _Nz;

  int64_t fieldSize{0};
  if (Nz > 1)
    fieldSize = (Nx + 2) * (Ny + 2) * (Nz + 2);
  else
    fieldSize = (Nx + 2) * (Ny + 2) * Nz;

  field.assign(fieldSize, 0.0);  // field.resize(Nx * Ny * Nz, 0.0);
}

Field::ComputingField::ComputingField(const ComputingField& _field) : field{_field.get_field()}
{
  field = _field.get_field();
  Nx = _field.get_Nx();
  Ny = _field.get_Ny();
  Nz = _field.get_Nz();
}

Field::ComputingField::ComputingField(ComputingField&& _field) noexcept : field{std::move(_field.get_field())}
{
  Nx = _field.get_Nx();
  Ny = _field.get_Ny();
  Nz = _field.get_Nz();
}

void Field::ComputingField::resize_field(const int64_t _Nx, const int64_t _Ny, const int64_t _Nz)
{
  // this->Field::ComputingField::ComputingField(_Nx, _Ny, _Nz);
  if (_Nx <= 0 || _Ny <= 0 || _Nz < 0)
  {
    //throw std::runtime_error("Error: The grid size is incorrect");
    std::cout << "\nError: The grid size is incorrect\n";
    exit(-1);
  }
  Nx = _Nx;
  Ny = _Ny;
  Nz = _Nz;
  field.assign((_Nx + 2) * (_Ny + 2) * _Nz, 0.0);
}

double& Field::ComputingField::operator()(const int64_t i, const int64_t j, const int64_t k) {

  return field[(Nx + 2) * (Ny + 2) * (k + 1) + (Nx + 2) * (j + 1) + (i + 1)];
  // return field[(Nx + 2) * (Ny + 2) * (k + 1) + (Ny + 2) * (i + 1) + (j + 1)];
}

const double& Field::ComputingField::operator()(const int64_t i, const int64_t j, const int64_t k) const
{
  return field[(Nx + 2) * (Ny + 2) * (k + 1) + (Nx + 2) * (j + 1) + (i + 1)];
  // return field[(Nx + 2) * (Ny + 2) * (k + 1) + (Ny + 2) * (i + 1) + (j + 1)];
}

Field::ComputingField& Field::ComputingField::operator=(const Field::ComputingField& _field) {
  if (&_field != this) {
    // this->Field::ComputingField::ComputingField(_field);
    Nx = _field.get_Nx();
    Ny = _field.get_Ny();
    Nz = _field.get_Nz();
    field = _field.get_field();
  }
  return *this;
}

Field::ComputingField& Field::ComputingField::operator=(ComputingField&& _field) noexcept
{
  if (this != &_field)
  {
    // this->Field::ComputingField::ComputingField(std::move(_field));
    Nx = _field.get_Nx();
    Ny = _field.get_Ny();
    Nz = _field.get_Nz();
    field = std::move(_field.get_field());
  }
  return *this;
}

// void Field::ComputingField::clear_file(const char* path)
// {
//   std::ofstream outfile;
//   outfile.open(path, std::ios::out | std::ios::trunc); // Очистка файла при открытии на запись
//   outfile.close();
// }

void Field::ComputingField::clear_files(const char *directory_path) {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  fs::path dir_path = fs::path(directory_path);

  // Формирование имени файла для текущего процесса
  // Обращение сначала к patent_path из-за того, что путь заканчивается на / (поэтому последний элемент пустой)
  fs::path file_path = fs::path(directory_path);
  if (file_path.parent_path().filename() == "my_data")
    file_path /= "my_data_" + std::to_string(rank) + ".csv";
  else if (file_path.parent_path().filename() == "analytical_data")
    file_path /= "analytical_data_" + std::to_string(rank) + ".csv";
  else {
    std::cerr << "Field::ComputingField::clear_files() error: invalid directory path" << std::endl;
    exit(-1);
  }

  // Очистка файла или его создание
  try {
      std::ofstream file(file_path, std::ofstream::out | std::ofstream::trunc);
      if (file.is_open()) {
          file.close();
      } else {
          std::cerr << "Process " << rank << " failed to open file: " << file_path << std::endl;
      }
  } catch (const std::exception& e) {
      std::cerr << "Error in process " << rank << " while handling file " << file_path << ": " << e.what() << std::endl;
  }

  // Удаление лишних файлов
  // В случае ручного вмешательства в директории, следует очистить файлы вручную
  if (rank == 0) {
    int i = size;
    while (true) {
      fs::path extra_file_path;
      if (dir_path.parent_path().filename() == "my_data")
        extra_file_path = dir_path / ("my_data_" + std::to_string(i) + ".csv");
      else if (dir_path.parent_path().filename() == "analytical_data")
        extra_file_path = dir_path / ("analytical_data_" + std::to_string(i) + ".csv");

      if (fs::exists(extra_file_path)) {
        try {
          fs::remove(extra_file_path);
        } catch (const std::exception &e) {
          std::cerr << "Error while removing extra file " << extra_file_path << ": " << e.what() << std::endl;
        }
      } else {
        // Прекращение цикла, если файл не существует
        // Так как файлы идут по порядку, то если один не существует, то и последующие тоже не существуют
        break; 
      }
      ++i;
    }
  }
}

// "index" - это индекс строки или столбца, который фиксируется.
void Field::ComputingField::write_field_to_file(const char* path, const int64_t index, Axis axis)
{
  std::ofstream outfile;
  outfile.open(path, std::ios::app);
  if (!outfile.is_open())
  {
    //throw std::runtime_error("The file can't be opened!");
    std::cout << "\nThe file can't be opened!" << std::endl;
    exit(-1);
  }

  auto check_index = [](int64_t _index, int64_t field_length) -> void
  {
    if (_index >= field_length)
    {
      //throw std::runtime_error("Error: Going beyond the column indexing in the matrix (Writing field to the file)");
      std::cout << "\nError: Going beyond the column indexing in the matrix (Writing field to the file)\n";
      exit(-1);
    }
  };

  switch (axis)
  {
   // Здесь 'index' пока для 2-мерного случая
  case Axis::Ox:
    check_index(index, Ny);
    for (int64_t i = 0; i < Nx; ++i)
      // outfile << this->operator()(i, index, -1) << (i < Nx - 1 ? ';' : '\n'); // Меняю с (i, index, 0) на (i, index, -1)
      outfile << this->operator()(i, index, 1) << (i < Nx - 1 ? ';' : '\n'); // Меняю с (i, index, 0) на (i, index, -1)
    break;
  case Axis::Oy:
    check_index(index, Nx);
    for (int64_t j = 0; j < Ny; ++j)
      outfile << this->operator()(index, j, -1) << (j < Ny - 1 ? ';' : '\n');
    break;
  case Axis::Oz:
    if (Nz == 1) return;
    check_index(index, Ny);
    for (int64_t k = 0; k < Nz; ++k) // k = 0, но по оси Oz двигаемся только в 3D, так что ОК
      outfile << this->operator()(0, index, k) << (k < Nz - 1 ? ';' : '\n');
    break;
  default:
    std::cout << "\nInvalid axis specified\n";
    exit(-1);
  }
  outfile.close();
}

void Field::ComputingField::write_field_to_file_OX(const char* path, const int64_t j)
{
  write_field_to_file(path, j, Axis::Ox);
}

void Field::ComputingField::write_field_to_file_OY(const char* path, const int64_t i)
{
  write_field_to_file(path, i, Axis::Oy);
}

void Field::ComputingField::write_field_to_file_OZ(const char* path, const int64_t j)
{
  write_field_to_file(path, j, Axis::Oz);
}

void Field::ComputingField::clear_field() noexcept {
  field.clear();
  Nx = Ny = Nz = 0;
}
