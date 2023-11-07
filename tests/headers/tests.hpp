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

  static double get_global_err(Field::ComputingField& field_1, Field::ComputingField& field_2);


  static void Courant_condition_check(const double _dt, const double _dx);

  TEST(Test, main_test) {
    std::pair<uint64_t, uint64_t> Nx_Ny = { 1000ull, 1000ull };
    std::pair<double, double> ax_ay = { 0.0, 0.0 };
    std::pair<double, double> bx_by = { 1.0, 1.0 };
    double dt = 2e-15;
    //double t = 2e-14;
    double t = 1e-14;

    FDTD::FDTD field(Nx_Ny, ax_ay, bx_by, dt);
    Courant_condition_check(dt, field.get_dx());

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

    Field::ComputingField::clear_file(path_to_calculated_data);
    Field::ComputingField::clear_file(path_to_analytic_data);

    field.write_fields_to_file(path_to_calculated_data);
    analytical_field.write_fields_to_file(path_to_analytic_data);
  }

  TEST(Test, Checking_the_convergence)
  {
    std::pair<uint64_t, uint64_t> Nx_Ny = { 1000ull, 1000ull };
    std::pair<double, double> ax_ay = { 0.0, 0.0 };
    std::pair<double, double> bx_by = { 1.0, 1.0 };

    double dt = 2e-15;
    //double t = 2e-14;
    double t = 1e-14;

    FDTD::FDTD field(Nx_Ny, ax_ay, bx_by, dt);
    Courant_condition_check(dt, field.get_dx());

    set_default_field(field, std::make_pair(ax_ay.first, bx_by.first));
    FDTD::FDTD analytical_field(field);

    analytical_default_solution(analytical_field,
      std::make_pair(ax_ay.first, bx_by.first), t);
    field.field_update(t);

    double first_err = get_global_err(analytical_field.get_Ez(), field.get_Ez());
    std::cout << "First_err = " << first_err << '\n';

    //=====Second field=====
    
    Nx_Ny = { 2000ull, 2000ull };

    FDTD::FDTD field_2(Nx_Ny, ax_ay, bx_by, dt);
    Courant_condition_check(dt, field_2.get_dx());

    set_default_field(field_2, std::make_pair(ax_ay.first, bx_by.first));
    FDTD::FDTD analytical_field_2(field_2);

    analytical_default_solution(analytical_field_2,
      std::make_pair(ax_ay.first, bx_by.first), t);
    field_2.field_update(t);

    double second_err = get_global_err(analytical_field_2.get_Ez(), field_2.get_Ez());
    std::cout << "Second_err = " << second_err << '\n';

  }
}  // namespace gtest

#endif  // !__TESTS_HPP__
