#ifndef __TESTS_HPP__
#define __TESTS_HPP__

#include "FDTD.hpp"
#include "gtest.h"

namespace gtest {

  static void set_default_field(FDTD::FDTD& _field,
    std::pair<double, double>& ax_bx);

  static void analytical_default_solution(FDTD::FDTD& _field,
    std::pair<double, double>& ax_bx,
    double t);

  TEST(Test, main_test) {
    std::pair<uint64_t, uint64_t> Nx_Ny = { 500ull, 510ull }; // { 1000ull, 1010ull }
    std::pair<double, double> ax_ay = { 0.0, 0.0 };
    std::pair<double, double> bx_by = { 1.0, 1.0 };
    double dt = 0.01;
    double t = 0.6;

    FDTD::FDTD field(Nx_Ny, ax_ay, bx_by, dt);
    set_default_field(field, std::make_pair(ax_ay.first, bx_by.first));
    FDTD::FDTD analytical_field(field);

    analytical_default_solution(analytical_field,
      std::make_pair(ax_ay.first, bx_by.first), t);
    field.field_update(t);

    for (uint64_t y = 0ull; y < field.get_Ny(); ++y)
      for (uint64_t x = 0ull; x < field.get_Nx(); ++x) {
        ASSERT_NEAR(analytical_field.get_Ex()(x, y), field.get_Ex()(x, y), 0.01);
        ASSERT_NEAR(analytical_field.get_Ey()(x, y), field.get_Ey()(x, y), 0.01);
        ASSERT_NEAR(analytical_field.get_Ez()(x, y), field.get_Ez()(x, y), 0.01);
        ASSERT_NEAR(analytical_field.get_Bx()(x, y), field.get_Bx()(x, y), 0.01);
        ASSERT_NEAR(analytical_field.get_By()(x, y), field.get_By()(x, y), 0.01);
        ASSERT_NEAR(analytical_field.get_Bz()(x, y), field.get_Bz()(x, y), 0.01);
      }
    field.get_Ey().write_to_file();
    analytical_field.get_Ey().write_to_file();
  }
}  // namespace gtest

#endif  // !__TESTS_HPP__
