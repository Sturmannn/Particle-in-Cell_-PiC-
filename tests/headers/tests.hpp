#ifndef __TESTS_HPP__
#define __TESTS_HPP__

#include "FDTD.hpp"
#include "gtest.h"

namespace gtest {

enum class Component { Ex, Ey, Ez, Bx, By, Bz };
enum class Shift { shifted, unshifted };

class Test_obj {
public:
  // Конструкторы и деструктор
  Test_obj() = delete;
  Test_obj(const Component _E, const Component _B, FDTD::FDTD &_field);
  Test_obj(const Component _E, const Component _B, FDTD::FDTD &&_field);
  Test_obj(const Test_obj &other_test_field);
  ~Test_obj() = default;

  // Операторы присваивания
  Test_obj &operator=(const Test_obj &other_test_field);
  Test_obj &operator=(const FDTD::FDTD &other_test_field);
  Test_obj &operator=(Test_obj &&other_test_field) noexcept;

  // Методы для аналитического решения и установки полей
  void analytical_default_solution_OX(const Component E, const Component B,
                                      const double t, const Shift _shift);
  void analytical_default_solution_OY(const Component E, const Component B,
                                      const double t, const Shift _shift);
  void set_default_field_OX(const Component E, const Component B,
                            const Shift _shift);
  void set_default_field_OY(const Component E, const Component B,
                            const Shift _shift);

  // Метод для дальнейшего вызова соответствующего численного решения (shifted /
  // unshifted)
  void numerical_solution(const double t, const Shift _shift);
  void numerical_solution(const uint64_t t, const Shift _shift);


  // Метод для вычисления глобальной ошибки
  double get_global_err(const Component component);

  // Метод для вывода результатов сходимости
  void print_convergence(Test_obj &other_test); // Check by component "E"
  // double get_convergence(); // Get by component "E"

  // Метод для записи сходимости в файл
  static void write_convergence_to_file(const char *path,
                                        std::vector<double> &data);

  FDTD::FDTD field;
  FDTD::FDTD analytical_field;

private:
  // Приватные члены данных
  Component E, B; // Компоненты поля

  // Проверка условия Куранта
  void Courant_condition_check(const Shift _shift) const noexcept;

  // Вспомогательные методы
  double helper_get_global_err(const Field::ComputingField &field_1,
                               const Field::ComputingField &field_2);
  double set_sign(const Component E, const Component B);
};

TEST(Test_version_comparison, shifted_OY) {
  std::pair<uint64_t, uint64_t> Nx_Ny = {64ull, 64ull};
  std::pair<double, double> ax_ay = {0.0, 0.0};
  std::pair<double, double> bx_by = {1.0, 1.0};
  double dt = 2e-15;
  // dt = 0.4625e-12;
  // double t = 2e-14;

  //double t = 1e-12;
  // t = 1e-10;

  double dx = (bx_by.first - ax_ay.first) / Nx_Ny.first;
  dt = 0.25 * dx / C;
  uint64_t t = 1000;

  Component E = Component::Ez;
  Component B = Component::Bx;
  Shift shift = Shift::unshifted;

  FDTD::FDTD field(Nx_Ny, ax_ay, bx_by, dt);
  Test_obj test(E, B, std::move(field));
  test.analytical_default_solution_OY(E, B, t * dt, shift);
  test.set_default_field_OY(E, B, shift);

  // test.analytical_default_solution_OX(E, B, t);
  // test.set_default_field_OX(E, B);
  test.numerical_solution(t, shift);

  Field::ComputingField::clear_file(path_to_calculated_data);
  Field::ComputingField::clear_file(path_to_analytic_data);

  test.field.write_fields_to_file_OY(path_to_calculated_data,
                                     test.field.get_dy());
  test.analytical_field.write_fields_to_file_OY(path_to_analytic_data,
                                                test.field.get_dy());

  // test.field.write_fields_to_file_OX(path_to_calculated_data,
  // test.field.get_dx());
  // test.analytical_field.write_fields_to_file_OX(path_to_analytic_data,
  // test.field.get_dx());
}

// TEST(Test, main_test) {
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

//      ASSERT_NEAR(analytical_field.get_Ex()(x, y), field.get_Ex()(x, y),
//      0.01); ASSERT_NEAR(analytical_field.get_Ey()(x, y), field.get_Ey()(x,
//      y), 0.01); ASSERT_NEAR(analytical_field.get_Ez()(x, y),
//      field.get_Ez()(x, y), 0.01); ASSERT_NEAR(analytical_field.get_Bx()(x,
//      y), field.get_Bx()(x, y), 0.01);
//      ASSERT_NEAR(analytical_field.get_By()(x, y), field.get_By()(x, y),
//      0.01); ASSERT_NEAR(analytical_field.get_Bz()(x, y), field.get_Bz()(x,
//      y), 0.01);
//    }

//  Field::ComputingField::clear_file(path_to_calculated_data);
//  Field::ComputingField::clear_file(path_to_analytic_data);

//  field.write_fields_to_file(path_to_calculated_data);
//  analytical_field.write_fields_to_file(path_to_analytic_data);
//}

// TEST(Test_version_comparison,
// UNSHIFTED_Checking_the_convergence__1_iteration)
//{
//  std::pair<uint64_t, uint64_t> Nx_Ny = { 64ull, 64ull };
//  std::pair<double, double> ax_ay = { 0.0, 0.0 };
//  std::pair<double, double> bx_by = { 1.0, 1.0 };
//
//  double dt = 1e-15;
//  //double t = 2e-14;
//  double t = 1e-13;
//
//  FDTD::FDTD field(Nx_Ny, ax_ay, bx_by, dt);
//
//
//  Component E = Component::Ez;
//  Component B = Component::Bx;
//  Shift shift = Shift::unshifted;
//
//  Test_obj test(E, B, std::move(field));
//  test.analytical_default_solution_OY(E, B, t, shift);
//  test.set_default_field_OY(E, B, shift);
//  test.numerical_solution(t, shift);
//
//
//
//  //Courant_condition_check(dt, field.get_dx());
//
//  //set_default_field(field, std::make_pair(ax_ay.first, bx_by.first));
//  //FDTD::FDTD analytical_field(field);
//
//  //analytical_default_solution(analytical_field,
//  //  std::make_pair(ax_ay.first, bx_by.first), t);
//  //field.field_update(t);
//
//  //double first_err = get_global_err(analytical_field.get_Ez(),
//  field.get_Ez());
//  //std::cout << "First_err = " << first_err << '\n';
//
//  //=====Second field=====
//
//  Nx_Ny = { 128ull, 128ull };
//  dt = dt / 4;
//
//  FDTD::FDTD field_2(Nx_Ny, ax_ay, bx_by, dt);
//
//
//  Test_obj test_2(E, B, std::move(field_2));
//  test_2.analytical_default_solution_OY(E, B, t, shift);
//  test_2.set_default_field_OY(E, B, shift);
//  test_2.numerical_solution(t, shift);
//
//  test.print_convergence(test_2);
//
//  //Courant_condition_check(dt, field_2.get_dx());
//
//  //set_default_field(field_2, std::make_pair(ax_ay.first, bx_by.first));
//  //FDTD::FDTD analytical_field_2(field_2);
//
//  //analytical_default_solution(analytical_field_2,
//  //  std::make_pair(ax_ay.first, bx_by.first), t);
//  //field_2.field_update(t);
//
//  //double second_err = get_global_err(analytical_field_2.get_Ez(),
//  field_2.get_Ez());
//  //std::cout << "Second_err = " << second_err << '\n';
//
//  //std::cout << "difference = " << first_err / second_err << '\n';
//
//}

TEST(Test_version_comparison, SHIFTED_Checking_the_convergence__1_iteration) {
  std::pair<uint64_t, uint64_t> Nx_Ny = {64ull, 64ull};
  std::pair<double, double> ax_ay = {0.0, 0.0};
  std::pair<double, double> bx_by = {1.0, 1.0};
  double start = omp_get_wtime();

  double dt = 1e-15;
  // double t = 2e-14;
  double dx = (bx_by.first - ax_ay.first) / Nx_Ny.first;
  dt = 0.25 * dx / C;
  //double t = 1e-12;


  //t = (bx_by.first - ax_ay.first) / C * 0.25;
  uint64_t t = 100; // Задание количества итераций

  FDTD::FDTD field(Nx_Ny, ax_ay, bx_by, dt);

  Component E = Component::Ez;
  Component B = Component::Bx;
  Shift shift = Shift::unshifted;

  Test_obj test(E, B, std::move(field));
  test.analytical_default_solution_OY(E, B, t * dt, shift);  // t * dt OR t
  test.set_default_field_OY(E, B, shift);
  test.numerical_solution(t, shift);

  //=====Second field=====

  Nx_Ny = {128ull, 128ull};
  dt = dt / 4;
  t *= 4;

  //dx = (bx_by.first - ax_ay.first) / Nx_Ny.first;
  //dt = 0.25 * dx / C;

  FDTD::FDTD field_2(Nx_Ny, ax_ay, bx_by, dt);

  Test_obj test_2(E, B, std::move(field_2));
  test_2.analytical_default_solution_OY(E, B, t * dt, shift); // t * dt OR t
  test_2.set_default_field_OY(E, B, shift);
  test_2.numerical_solution(t, shift);

  double end = omp_get_wtime();
  std::cout << "Time = " << end - start << std::endl;

  test.print_convergence(test_2);
}

// TEST(Test_version_comparison,
// SHIFTED_Checking_the_convergence__several_iterations)
//{
//  std::pair<uint64_t, uint64_t> Nx_Ny = { 64ull, 64ull };
//  std::pair<double, double> ax_ay = { 0.0, 0.0 };
//  std::pair<double, double> bx_by = { 1.0, 1.0 };

//  double dt = 1e-15;
//  //double t = 2e-14;
//  double t = 1e-13;

//  Component E = Component::Ez;
//  Component B = Component::Bx;
//  Shift shift = Shift::shifted;
//  double error = 0.0;

//  std::vector<double> convergences;
//  //FDTD::FDTD field(Nx_Ny, ax_ay, bx_by, dt);
//  //Test_obj test(E, B, std::move(field));
//  //test.analytical_default_solution_OY(E, B, t, shift);
//  //test.set_default_field_OY(E, B, shift);
//  //test.numerical_solution(t, shift);
//  //convergences.push_back(test.get_global_err(E));
//  for (uint64_t i = 0ull; i < 3ull; ++i)
//  {
//    Nx_Ny.first *= (1ull << i);
//    Nx_Ny.second *= (1ull << i);
//    dt /= pow(2, i);
//    //FDTD::FDTD field(Nx_Ny, ax_ay, bx_by, dt / (1ull >> i));
//    FDTD::FDTD field(Nx_Ny, ax_ay, bx_by, dt);
//    Test_obj test(E, B, field);
//    test.analytical_default_solution_OY(E, B, t, shift);
//    test.set_default_field_OY(E, B, shift);
//    test.numerical_solution(t, shift);

//    error = test.get_global_err(E);
//    convergences.push_back(error); // E
//    std::cout << error << ' ';
//    field.~FDTD();
//    test.~Test_obj();
//  }
//  std::cout << std::endl;
//  Test_obj::write_convergence_to_file(path_to_convergence_data, convergences);
//}
} // namespace gtest

#endif // !__TESTS_HPP__
