#ifndef __TESTS_HPP__
#define __TESTS_HPP__

#include "FDTD.hpp"
#include "gtest.h"

namespace gtest {

  enum class Component {Ex, Ey, Ez, Bx, By, Bz};

  class Test_obj
  {
  public:
    Test_obj() = delete;
    Test_obj(const Component _E, const Component _B, FDTD::FDTD&& _field);
    Test_obj(const Test_obj& other_test_field);
    ~Test_obj() = default;

    Test_obj& operator = (const Test_obj& other_test_field);
    Test_obj& operator = (Test_obj&& other_test_field) noexcept;

    void analytical_default_solution_OX(const Component E, const Component B, const double t);
    void analytical_default_solution_OY(const Component E, const Component B, const double t);
    void set_default_field_OX(const Component E, const Component B);
    void set_default_field_OY(const Component E, const Component B);
    void numerical_solution(const double t);
    double get_global_err(const Component component);
    void check_ñonvergence(Test_obj& other_test); // Check by component "E"

    FDTD::FDTD field;
    FDTD::FDTD analytical_field;
  private:
    Component E, B;
    void Courant_condition_check() const noexcept;
    double helper_get_global_err(const Field::ComputingField& _E, const Field::ComputingField& _B);
    double set_sign(const Component E, const Component B);
  };

  static void set_default_field_OY(FDTD::FDTD& _field, const Component E, const Component B);

  static void set_default_field(FDTD::FDTD& _field,
    std::pair<double, double>& ax_bx);

  static void analytical_default_solution(FDTD::FDTD& _field,
    std::pair<double, double>& ax_bx,
    double t);

  static double get_global_err(Field::ComputingField& field_1, Field::ComputingField& field_2);


  static void Courant_condition_check(const double _dt, const double _dx);

  TEST(Test_construcrot, proverka)
  {
    std::pair<uint64_t, uint64_t> Nx_Ny = { 64ull, 64ull };
    std::pair<double, double> ax_ay = { 0.0, 0.0 };
    std::pair<double, double> bx_by = { 1.0, 1.0 };
    double dt = 2e-15;
    //double t = 2e-14;
    double t = 1e-12;

    Component E = Component::Ez;
    Component B = Component::Bx;

    FDTD::FDTD field(Nx_Ny, ax_ay, bx_by, dt);
    Test_obj test(E, B, std::move(field));
    test.analytical_default_solution_OY(E, B, t);
    test.set_default_field_OY(E, B);

    //test.analytical_default_solution_OX(E, B, t);
    //test.set_default_field_OX(E, B);
    test.numerical_solution(t);

    Field::ComputingField::clear_file(path_to_calculated_data);
    Field::ComputingField::clear_file(path_to_analytic_data);

    test.field.write_fields_to_file_OY(path_to_calculated_data, test.field.get_dy());
    test.analytical_field.write_fields_to_file_OY(path_to_analytic_data, test.field.get_dy());

    //test.field.write_fields_to_file_OX(path_to_calculated_data, test.field.get_dx());
    //test.analytical_field.write_fields_to_file_OX(path_to_analytic_data, test.field.get_dx());
  }

  //TEST(Test, main_test) {
  //  std::pair<uint64_t, uint64_t> Nx_Ny = { 64ull, 64ull };
  //  std::pair<double, double> ax_ay = { 0.0, 0.0 };
  //  std::pair<double, double> bx_by = { 1.0, 1.0 };
  //  double dt = 2e-15;
  //  //double t = 2e-14;
  //  double t = 1e-12;

  //  FDTD::FDTD field(Nx_Ny, ax_ay, bx_by, dt);
  //  Courant_condition_check(dt, field.get_dx());

  //  set_default_field(field, std::make_pair(ax_ay.first, bx_by.first));
  //  FDTD::FDTD analytical_field(field);

  //  analytical_default_solution(analytical_field,
  //    std::make_pair(ax_ay.first, bx_by.first), t);
  //  field.field_update(t);

  //  for (uint64_t y = 0ull; y < field.get_Ny(); ++y)
  //    for (uint64_t x = 0ull; x < field.get_Nx(); ++x) {

  //      ASSERT_NEAR(analytical_field.get_Ex()(x, y), field.get_Ex()(x, y), 0.01);
  //      ASSERT_NEAR(analytical_field.get_Ey()(x, y), field.get_Ey()(x, y), 0.01);
  //      ASSERT_NEAR(analytical_field.get_Ez()(x, y), field.get_Ez()(x, y), 0.01);
  //      ASSERT_NEAR(analytical_field.get_Bx()(x, y), field.get_Bx()(x, y), 0.01);
  //      ASSERT_NEAR(analytical_field.get_By()(x, y), field.get_By()(x, y), 0.01);
  //      ASSERT_NEAR(analytical_field.get_Bz()(x, y), field.get_Bz()(x, y), 0.01);
  //    }

  //  Field::ComputingField::clear_file(path_to_calculated_data);
  //  Field::ComputingField::clear_file(path_to_analytic_data);

  //  field.write_fields_to_file(path_to_calculated_data);
  //  analytical_field.write_fields_to_file(path_to_analytic_data);
  //}

  TEST(Test, Checking_the_convergence)
  {
    std::pair<uint64_t, uint64_t> Nx_Ny = { 64ull, 64ull };
    std::pair<double, double> ax_ay = { 0.0, 0.0 };
    std::pair<double, double> bx_by = { 1.0, 1.0 };

    double dt = 1e-15;
    //double t = 2e-14;
    double t = 1e-13;

    FDTD::FDTD field(Nx_Ny, ax_ay, bx_by, dt);


    Component E = Component::Ez;
    Component B = Component::Bx;

    Test_obj test(E, B, std::move(field));
    test.analytical_default_solution_OY(E, B, t);
    test.set_default_field_OY(E, B);
    test.numerical_solution(t);



    //Courant_condition_check(dt, field.get_dx());

    //set_default_field(field, std::make_pair(ax_ay.first, bx_by.first));
    //FDTD::FDTD analytical_field(field);

    //analytical_default_solution(analytical_field,
    //  std::make_pair(ax_ay.first, bx_by.first), t);
    //field.field_update(t);

    //double first_err = get_global_err(analytical_field.get_Ez(), field.get_Ez());
    //std::cout << "First_err = " << first_err << '\n';

    //=====Second field=====
    
    Nx_Ny = { 128ull, 128ull };
    dt = dt / 4;

    FDTD::FDTD field_2(Nx_Ny, ax_ay, bx_by, dt);


    Test_obj test_2(E, B, std::move(field_2));
    test_2.analytical_default_solution_OY(E, B, t);
    test_2.set_default_field_OY(E, B);
    test_2.numerical_solution(t);

    test.check_ñonvergence(test_2);

    //Courant_condition_check(dt, field_2.get_dx());

    //set_default_field(field_2, std::make_pair(ax_ay.first, bx_by.first));
    //FDTD::FDTD analytical_field_2(field_2);

    //analytical_default_solution(analytical_field_2,
    //  std::make_pair(ax_ay.first, bx_by.first), t);
    //field_2.field_update(t);

    //double second_err = get_global_err(analytical_field_2.get_Ez(), field_2.get_Ez());
    //std::cout << "Second_err = " << second_err << '\n';

    //std::cout << "difference = " << first_err / second_err << '\n';

  }
}  // namespace gtest

#endif  // !__TESTS_HPP__
