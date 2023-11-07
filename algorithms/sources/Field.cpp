#include "Field.hpp"

Field::ComputingField::ComputingField(const uint64_t _Nx, const uint64_t _Ny)
  : Nx(_Nx), Ny(_Ny) {
  field.resize(Nx * Ny, 0.0);
}

Field::ComputingField::ComputingField(const ComputingField& _field)
  : field(_field.field), Nx(_field.Nx), Ny(_field.Ny) {}

double& Field::ComputingField::operator()(uint64_t i, uint64_t j) {

  //if (i == SIZE_MAX) i = 0ull;
  //if (j == SIZE_MAX) j = 0ull;

  if (i == SIZE_MAX) i = Nx - 1ull;
  if (j == SIZE_MAX) j = Ny - 1ull;

  return field[Nx * (j % Ny) + (i % Nx)];
}

//double& Field::ComputingField::operator()(const int64_t i, const int64_t j) {
//  
//  int64_t x, y;
//  if (i >= Nx) x = 0;
//  else if (i < 0) x = Nx - 1;
//  else x = i;
//
//  if (i >= Ny) y = 0;
//  else if (i < 0) y = Ny - 1;
//  else y = i;
//
//  return field[Nx * y + x];
//  //return field[Nx * (j % Ny) + (i % Nx)];
//}

Field::ComputingField& Field::ComputingField::operator=(
  const Field::ComputingField& _field) {
  if (&_field != this) {
    Nx = _field.Nx;
    Ny = _field.Ny;
    field.resize(Nx * Ny);
    field.assign(_field.field.begin(), _field.field.end());
  }
  return *this;
}

void Field::ComputingField::clear_file(const char* path)
{
  std::ofstream outfile;
  outfile.open(path, std::ios::out | std::ios::trunc);
  outfile.close();
}

void Field::ComputingField::write_field_to_file(const char* path, const uint64_t j)
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
