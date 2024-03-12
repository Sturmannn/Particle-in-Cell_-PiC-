#include "Field.hpp"

Field::ComputingField::ComputingField(const uint64_t _Nx, const uint64_t _Ny, const uint64_t _Nz)
{
  Nx = _Nx;
  Ny = _Ny;

  if (_Nz == 0)
  {
    Nz = 1;
    field.resize(Nx * Ny, 0.0);
  }
  else
  {
    Nz = _Nz;
    field.resize(Nx * Ny * Nz, 0.0);
  }
}

Field::ComputingField::ComputingField(const ComputingField& _field)
{
  field = _field.field;
  Nx = _field.get_Nx();
  Ny = _field.get_Ny();
  Nz = _field.get_Nz();
}

Field::ComputingField::ComputingField(ComputingField&& _field) noexcept : field(std::move(_field.field))
{
  Nx = _field.get_Nx();
  Ny = _field.get_Ny();
  Nz = _field.get_Nz();
}

void Field::ComputingField::resize_field(const uint64_t _Nx, const uint64_t _Ny, const uint64_t _Nz)
{
  this->Field::ComputingField::ComputingField(_Nx, _Ny, _Nz);
}

double& Field::ComputingField::operator()(uint64_t i, uint64_t j, uint64_t k) {

  //if (i == SIZE_MAX) i = 0ull;
  //if (j == SIZE_MAX) j = 0ull;

  if (i == SIZE_MAX) i = Nx - 1ull;
  if (j == SIZE_MAX) j = Ny - 1ull;
  if (k == SIZE_MAX) k = Nz - 1ull;

  //return field[Nx * (j % Ny) + (i % Nx)];
  return field[Nx * Ny * (k % Nz) + Nx * (j % Ny) + (i % Nx)];
}

const double& Field::ComputingField::operator()(uint64_t i, uint64_t j, uint64_t k) const
{
  //return const_cast<const double&>((*this)(i, j, k));
  //return (*this)(i, j, k);

  if (i == SIZE_MAX) i = Nx - 1ull;
  if (j == SIZE_MAX) j = Ny - 1ull;
  if (k == SIZE_MAX) k = Nz - 1ull;

  return field[Nx * Ny * (k % Nz) + Nx * (j % Ny) + (i % Nx)];
}

Field::ComputingField& Field::ComputingField::operator=(
  const Field::ComputingField& _field) {
  if (&_field != this) {
    this->Field::ComputingField::ComputingField(_field);
  }
  return *this;
}

Field::ComputingField& Field::ComputingField::operator=(ComputingField&& _field) noexcept
{
  if (this != &_field)
  {
    this->Field::ComputingField::ComputingField(std::move(_field));
  }
  return *this;
}

void Field::ComputingField::clear_file(const char* path)
{
  std::ofstream outfile;
  outfile.open(path, std::ios::out | std::ios::trunc);
  outfile.close();
}

void Field::ComputingField::write_field_to_file_OX(const char* path, const uint64_t j)
{
  std::ofstream outfile;
  outfile.open(path, std::ios::app);
  if(!outfile.is_open())
  {
    std::cout << "The file can't be opened!" << std::endl;
    exit(-1);
  }
  if (j >= Ny)
  {
    std::cout << "Error: Going beyond the column indexing in the matrix (Writing field to the file)\n";
    exit(-1);
  }
  uint64_t i = 0ull;
  for (; i < Nx - 1ull; ++i)
    outfile << this->operator()(i, j) << ';';
  outfile << this->operator()(i, j) << '\n';
  outfile.close();
}

void Field::ComputingField::write_field_to_file_OY(const char* path, const uint64_t i)
{
  std::ofstream outfile;
  outfile.open(path, std::ios::app);
  if (!outfile.is_open())
  {
    std::cout << "The file can't be opened!" << std::endl;
    exit(-1);
  }
  if (i >= Nx)
  {
    std::cout << "Error: Going beyond the column indexing in the matrix (Writing field to the file)\n";
    exit(-1);
  }
  uint64_t j = 0ull;
  for (; j < Ny - 1ull; ++j)
    outfile << this->operator()(i, j) << ';';
  outfile << this->operator()(i, j) << '\n';
  outfile.close();
}

void Field::ComputingField::write_field_to_file_OZ(const char* path, uint64_t k)
{
  std::ofstream outfile;
  outfile.open(path, std::ios::app);
  if (!outfile.is_open())
  {
    std::cout << "The file can't be opened!" << std::endl;
    exit(-1);
  }
  if (k >= Nz)
  {
    std::cout << "Error: Going beyond the column indexing in the matrix (Writing field to the file)\n";
    exit(-1);
  }
  uint64_t j = 0ull;
  uint64_t i = 0ull;
  k = 0;
  for (; k < Nz - 1ull; ++k)
    outfile << this->operator()(i, j, k) << ';';
  outfile << this->operator()(i, j, k) << '\n';
  outfile.close();
}
