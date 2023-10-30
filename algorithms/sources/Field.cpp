#include "Field.hpp"

Field::ComputingField::ComputingField(const uint64_t _Nx, const uint64_t _Ny)
  : Nx(_Nx), Ny(_Ny) {
  field.resize(Nx * Ny, 0.0);
}

Field::ComputingField::ComputingField(const ComputingField& _field)
  : field(_field.field), Nx(_field.Nx), Ny(_field.Ny) {}

double& Field::ComputingField::operator()(const uint64_t i, const uint64_t j) {
  return field[Nx * (j % Ny) + (i % Nx)];
}

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

void Field::ComputingField::write_to_file(const uint64_t j)
{
  std::ofstream an_outfile(path_to_analytic_data);
  std::ofstream my_outfile(path_to_calculated_data);
  if (!an_outfile.is_open() || !my_outfile.is_open())
  {
    std::cout << "One of the files can't be opened!\n";
    exit(-1);
  }
  for (uint64_t i = 0ull; i < Nx; ++i)
  {
    an_outfile << this->operator()(i, j) << ';';
  
    my_outfile << this->operator()(i, j) << ';';
  }
  
  std::cout.flush();
  an_outfile.close();
  my_outfile.close();
}
