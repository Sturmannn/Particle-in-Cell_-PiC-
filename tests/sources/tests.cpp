#include "tests.hpp"

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

void gtest::set_default_field(FDTD::FDTD& _field,
  std::pair<double, double>& ax_bx) {
  Field::ComputingField& Ex = _field.get_Ex();
  Field::ComputingField& Ey = _field.get_Ey();
  Field::ComputingField& Ez = _field.get_Ez();

  Field::ComputingField& Bx = _field.get_Bx();
  Field::ComputingField& By = _field.get_By();
  Field::ComputingField& Bz = _field.get_Bz();

  for (uint64_t y = 0; y < _field.get_Ny(); ++y)
    for (uint64_t x = 0; x < _field.get_Nx(); ++x) {
      Ey(x, y) = Bz(x, y) =
        sin(2.0 * PI * (static_cast<double>(x) - ax_bx.first) /
          (ax_bx.second - ax_bx.first));
      Ex(x, y) = Ez(x, y) = Bx(x, y) = By(x, y) = 0.0;
    }
}

void gtest::analytical_default_solution(FDTD::FDTD& _field,
  std::pair<double, double>& ax_bx,
  double t) {
  Field::ComputingField& Ex = _field.get_Ex();
  Field::ComputingField& Ey = _field.get_Ey();
  Field::ComputingField& Ez = _field.get_Ez();

  Field::ComputingField& Bx = _field.get_Bx();
  Field::ComputingField& By = _field.get_By();
  Field::ComputingField& Bz = _field.get_Bz();

  for (uint64_t y = 0ull; y < _field.get_Ny(); ++y)
    for (uint64_t x = 0ull; x < _field.get_Nx(); ++x) {
      Ey(x, y) = Bz(x, y) =
        sin(2.0 * PI * (static_cast<double>(x) - ax_bx.first - C * t) /
          (ax_bx.second - ax_bx.first));
      Ex(x, y) = Ez(x, y) = Bx(x, y) = By(x, y) = 0.0;
    }
}
