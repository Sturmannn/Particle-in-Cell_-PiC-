#include "Field.hpp"

Field::ComputingField::ComputingField(const uint64_t _Nx, const uint64_t _Ny)
  : Nx(_Nx), Ny(_Ny) {
  field.resize(Nx * Ny, 0.0);
}

Field::ComputingField::ComputingField(const ComputingField& _field)
  : field(_field.field), Nx(_field.Nx), Ny(_field.Ny) {}

double& Field::ComputingField::operator()(uint64_t i, uint64_t j) {

  if (i == SIZE_MAX) i = 0ull;
  if (j == SIZE_MAX) j = 0ull;

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

void Field::ComputingField::write_to_file(const double _dx, const uint64_t j)
{
  std::ofstream an_outfile(path_to_analytic_data);
  std::ofstream my_outfile(path_to_calculated_data);
  if (!an_outfile.is_open() || !my_outfile.is_open())
  {
    std::cout << "One of the files can't be opened!" << std::endl;
    exit(-1);
  }
  
  uint64_t i = 0ull;
  for (; i < Nx - 1ull; ++i)
  {
    an_outfile << this->operator()(i, j) << ';';
  
    my_outfile << this->operator()(i, j) << ';';
  }

  an_outfile << this->operator()(i, j) << '\n' << _dx;

  my_outfile << this->operator()(i, j) << '\n' << _dx;
  
  an_outfile.close();
  my_outfile.close();
}
