#include "tests.hpp"

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

void gtest::set_default_field(FDTD::FDTD& _field,
  std::pair<double, double>& ax_bx) {

  double x = ax_bx.first;

  Field::ComputingField& Ex = _field.get_Ex();
  Field::ComputingField& Ey = _field.get_Ey();
  Field::ComputingField& Ez = _field.get_Ez();

  Field::ComputingField& Bx = _field.get_Bx();
  Field::ComputingField& By = _field.get_By();
  Field::ComputingField& Bz = _field.get_Bz();

  for (uint64_t i = 0; i < _field.get_Nx(); ++i, x += _field.get_dx())
    for (uint64_t j = 0; j < _field.get_Ny(); ++j) {
      Ey(i, j) = Bz(i, j) =
        sin(2.0 * PI * (x - ax_bx.first) /
          (ax_bx.second - ax_bx.first));
      Ex(i, j) = Ez(i, j) = Bx(i, j) = By(i, j) = 0.0;
    }
}

void gtest::analytical_default_solution(FDTD::FDTD& _field,
  std::pair<double, double>& ax_bx,
  double t) {

  double x = ax_bx.first;

  Field::ComputingField& Ex = _field.get_Ex();
  Field::ComputingField& Ey = _field.get_Ey();
  Field::ComputingField& Ez = _field.get_Ez();

  Field::ComputingField& Bx = _field.get_Bx();
  Field::ComputingField& By = _field.get_By();
  Field::ComputingField& Bz = _field.get_Bz();

  for (uint64_t i = 0ull; i < _field.get_Nx(); ++i, x += _field.get_dx())
    for (uint64_t j = 0ull; j < _field.get_Ny(); ++j) {
      Ey(i, j) = Bz(i, j) =
        sin(2.0 * PI * (x - ax_bx.first - C * t) /
          (ax_bx.second - ax_bx.first));
      Ex(i, j) = Ez(i, j) = Bx(i, j) = By(i, j) = 0.0;
    }
}
