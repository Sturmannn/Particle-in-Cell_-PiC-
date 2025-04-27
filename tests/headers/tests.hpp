#ifndef __TESTS_HPP__
#define __TESTS_HPP__

#include <functional>
#include <fstream>
#include "FDTD.hpp"
#include <gtest/gtest.h>
#include <omp.h>


//using FDTD::Axis;
using FDTD::Component;

constexpr const char *path_to_measurements_file = PATH_TO_MEASUREMENTS_FILE;

namespace gtest {

extern int arg_Nx, arg_Ny, arg_Nz;

enum class Shift { shifted, unshifted };

class Test_obj {
public:
  Test_obj() = delete;
  Test_obj(const Component _E, const Component _B, int64_t Nx, int64_t Ny, int64_t Nz, int MPI_dimension); // Вариант пока что для 2D
  Test_obj(const Component _E, const Component _B, FDTD::FDTD &_field);
  Test_obj(const Component _E, const Component _B, FDTD::FDTD &&_field);
  Test_obj(const Test_obj &other_test_field);
  ~Test_obj() {MPI_Comm_free(&cart_comm);};

  void initialize_field(const FDTD::FDTD &field);
  void initialize_field(FDTD::FDTD &&field);

  Test_obj &operator=(const Test_obj &other_test_field);
  Test_obj &operator=(const FDTD::FDTD &other_test_field);
  Test_obj &operator=(Test_obj &&other_test_field) noexcept;

  FDTD::FDTD &get_field() { return field; }
  FDTD::FDTD &get_analytical_field() { return analytical_field; }
  const FDTD::FDTD &get_field() const { return field; }
  const FDTD::FDTD &get_analytical_field() const { return analytical_field; }

  void analytical_default_solution(const Component E, const Component B,
                                   const double t, const Shift _shift);

  void set_default_field(const Component E, const Component B,
                         const Shift _shift);

  void numerical_solution(const double t, const Shift _shift, MPI_Comm cart_comm);
  void numerical_solution(const int64_t t, const Shift _shift, MPI_Comm cart_comm);

  double get_delta_space(void) const;

  double get_global_err(const Component component);

  void print_convergence(Test_obj &other_test); // Check by component "E"

  static void write_convergence_to_file(const char *path,
                                        std::vector<double> &data);

  MPI_Comm &get_cart_comm(void) { return cart_comm; }
  
  std::vector<int>& get_coords() { return coords; }
  std::vector<int>& get_dims() { return dims; }
  // int (&get_dims(void))[3] { return dims; }

  int64_t get_local_Nx(void) const noexcept { return local_Nx; }
  int64_t get_local_Ny(void) const noexcept { return local_Ny; }
  int64_t get_local_Nz(void) const noexcept { return local_Nz; }
  std::tuple<int64_t, int64_t, int64_t> get_local_Nx_Ny_Nz(void) const noexcept {
    return std::make_tuple(local_Nx, local_Ny, local_Nz);
  }



private:
  FDTD::FDTD field;
  FDTD::FDTD analytical_field;
  Component E, B; // Компоненты поля
  MPI_Comm cart_comm; // Декартова топология
  std::vector<int> dims, coords; // Размеры топологии и координаты процесса (или 2 или 3 элемента)
  int64_t local_Nx, local_Ny, local_Nz; // Размеры поддомена

  void Courant_condition_check(const Shift _shift) const noexcept;
  void set_default_field_or_analytical_default_solution(const Component E,
                                                        const Component B,
                                                        const Shift _shift,
                                                        std::function<void(std::tuple<Axis, int64_t, int64_t>,
                                                          std::tuple<Axis, int64_t, int64_t>,
                                                          std::tuple<Axis, int64_t, int64_t>)> loop_function);
  void create_cartesian_topology(int world_size, int MPI_dimension);
  
  // Пока что здесь только 2D случай
  void set_subdomain_sizes(int64_t Nx, int64_t Ny, int64_t Nz);
};

TEST(Test_version_comparison, shifted_OZ) {

  // При изменении размеров сетки, не забыть изменить и размеры поддомена
  int MPI_dimension = 3;
  std::tuple<int64_t, int64_t, int64_t> Nx_Ny_Nz = {gtest::arg_Nx, gtest::arg_Ny, gtest::arg_Nz};
  std::tuple<double, double, double> ax_ay_az = {0.0, 0.0, 0.0};
  std::tuple<double, double, double> bx_by_bz = {1.0, 1.0, 1.0};
  std::tuple<double, double, double> dx_dy_dz = 
        {(std::get<0>(bx_by_bz) - std::get<0>(ax_ay_az)) / static_cast<double>(std::get<0>(Nx_Ny_Nz)),
         (std::get<1>(bx_by_bz) - std::get<1>(ax_ay_az)) / static_cast<double>(std::get<1>(Nx_Ny_Nz)),
         (std::get<2>(bx_by_bz) - std::get<2>(ax_ay_az)) / static_cast<double>(std::get<2>(Nx_Ny_Nz))};
  
  // if (std::get<2>(Nx_Ny_Nz) > 1) {
  //   std::get<2>(dx_dy_dz) = 0.0;
  // }
  // double dt = 2e-15;
  // omp_set_num_threads(3);
  printf("Max threads num: %d\n", omp_get_max_threads());
  int thread_num = 0;
  #pragma omp parallel
  {
    thread_num = omp_get_thread_num();
    printf("Hello from thread %d\n", thread_num);
    if (omp_get_num_threads() == 0) {
      printf("Number of threads: %d\n", omp_get_num_threads());
    }
  }


  double dx = std::get<0>(dx_dy_dz);
  double dt = 0.25 * dx / C;
  // int64_t t = 250; // Задание количества итераций
  int64_t t = 150; // Задание количества итераций

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Component E = Component::Ez;
  Component B = Component::Bx;
  Shift shift = Shift::shifted;
  
  if (rank == 0) {
    std::cout << "Nx = " << std::get<0>(Nx_Ny_Nz) << std::endl << \
      "Ny = "<< std::get<1>(Nx_Ny_Nz) << std::endl << \
      "Nz = " << std::get<2>(Nx_Ny_Nz) << std::endl;
    std::cout << "Axis: " << FDTD::FDTD::axisToString(E, B) << std::endl;
  }


// =========================================================


  Test_obj test(E, B, std::get<0>(Nx_Ny_Nz), std::get<1>(Nx_Ny_Nz), std::get<2>(Nx_Ny_Nz), MPI_dimension);
  // Внимательно проверять, что в field нужно передавать именно local размеры
  FDTD::FDTD field(test.get_local_Nx_Ny_Nz(), ax_ay_az, bx_by_bz, dx_dy_dz, dt); 
  test.initialize_field(field);

  test.analytical_default_solution(E, B, t * dt, shift);
  test.set_default_field(E, B, shift);
  
  MPI_Barrier(MPI_COMM_WORLD);
  double start_time = MPI_Wtime();
  test.numerical_solution(t, shift, test.get_cart_comm());
  MPI_Barrier(MPI_COMM_WORLD);
  double end_time = MPI_Wtime();
  double elapsed_time = end_time - start_time;
  if (rank == 0) {
    std::cout << "Elapsed time: " << elapsed_time << " seconds" << std::endl;
    
    std::vector<int> dims = test.get_dims();
    std::ofstream measure_file;
    
    measure_file.open(path_to_measurements_file, std::ios::app);
    if (measure_file.is_open()) {
      // Здесь только для 3D для замеров времени
      measure_file << elapsed_time << '\t' << omp_get_max_threads() << '\t' << world_size << "\t" <<
        dims[0] << '\t' << dims[1] << '\t' << dims[2] << '\t' << gtest::arg_Nx  << '\t' <<
          gtest::arg_Ny << '\t' << gtest::arg_Nz << std::endl;
    }
    else {
      std::cerr << "The file for measurements can't be opened" << std::endl;
    }
    measure_file.close();
  } 

  Field::ComputingField::clear_files(path_to_calculated_data_directory);
  Field::ComputingField::clear_files(path_to_analytical_data_directory);
  

  // Поскольку MPI разбиение на данный момент только по OY, то в python скрипте для графиков происходит сбор всех данных \
  с соответствующих процессов. Поэтому, при запуске на OX, данные дублируются в графиках
  test.get_field().write_fields_to_file(path_to_calculated_data_directory, E, B,
                                  test.get_delta_space());
    std::cout << "Hello from rank " << rank << ", "<< test.get_coords()[0] << " " << test.get_coords()[1] << std::endl;
  test.get_analytical_field().write_fields_to_file(path_to_analytical_data_directory, E, B,
                                              test.get_delta_space());

// =============================================================================

  ////// Проверка сходимости, создаём новый объект:
  // НЕ ЗАБЫТЬ! Если уменьшяю dt в 4 раза, то и t увеличиваю в 4 раза
  Nx_Ny_Nz = {128, 128, 1};
  // Ошибка уменьшается в 4 раза
  if (shift == Shift::shifted) { 
    dt /= 2;
    t *= 2;
  }
  // Ошибка уменьшается в 4 раза 
  else if (shift == Shift::unshifted) {
    dt /= 4;
    t *= 4;
  }
    // double start_time = omp_get_wtime();
    // double end_time = omp_get_wtime();
  // std::cout << "Time = " << end_time - start_time << std::endl;


// =============================================================================

  // local_Nx = std::get<0>(Nx_Ny_Nz);
  // local_Ny = (std::get<1>(Nx_Ny_Nz) + 2 * world_size) / world_size;

  // // Это рассчитывается с учётом дополнительных граничных полей.
  // // В реализации же этого кода, они учитываются в конструктуре, поэтому этот результат нужно будет уменьшить на 2
  // if (rank < (std::get<1>(Nx_Ny_Nz) + 2 * world_size) % world_size) local_Ny += 1;
  
  // // Уменьшаем на 2, так как в конструкторе учитываются граничные поля
  // // if (world_size > 1) local_Ny -= 2;
  // local_Ny -= 2;

  // // Пока что Nz = 1
  // local_Nx_Ny_Nz = {local_Nx, local_Ny, 1};

  // dx_dy_dz = 
  //       {(std::get<0>(bx_by_bz) - std::get<0>(ax_ay_az)) / static_cast<double>(std::get<0>(Nx_Ny_Nz)),
  //        (std::get<1>(bx_by_bz) - std::get<1>(ax_ay_az)) / static_cast<double>(std::get<1>(Nx_Ny_Nz)),
  //        (std::get<2>(bx_by_bz) - std::get<2>(ax_ay_az)) / static_cast<double>(std::get<2>(Nx_Ny_Nz))};

  // FDTD::FDTD field_2(local_Nx_Ny_Nz, ax_ay_az, bx_by_bz, dx_dy_dz, dt);
  // Test_obj other_test(E,B, std::move(field_2));
  // other_test.analytical_default_solution(E, B, t * dt, shift);
  // other_test.set_default_field(E, B, shift);
  // other_test.numerical_solution(t, shift); // t для без сдвигов, а t * dt для сдвигов (не так...Учёт идёт в int/double)
  // test.print_convergence(other_test);

// =============================================================================
}







// TEST(Test, unshifted_field_test_several_iterations_2d) {
//   std::tuple<int64_t, int64_t, int64_t> Nx_Ny_Nz = {16, 16, 1}; // Start size of the field
//   std::tuple<double, double, double> ax_ay_az = {0.0, 0.0, 0.0};
//   std::tuple<double, double, double> bx_by_bz = {1.0, 1.0, 1.0};

//   auto increase_Nx_Ny_Nz = [](std::tuple<int64_t, int64_t, int64_t> &_Nx_Ny_Nz, int64_t value) -> void {
//     std::get<0>(_Nx_Ny_Nz) *= value;
//     std::get<1>(_Nx_Ny_Nz) *= value;
//     if (std::get<2>(_Nx_Ny_Nz) > 2) {
//       std::get<2>(_Nx_Ny_Nz) *= value;
//     }
//   };

//   // dx здесь для того, чтобы вычислить dt
//   double dx = (std::get<0>(bx_by_bz) - std::get<0>(ax_ay_az)) / std::get<0>(Nx_Ny_Nz);
//   double dt = 0.25 * dx / C;
//   int64_t t = 10; // Задание количества итераций

//   // Положительное направление ОУ
//   Component E = Component::Ex;
//   Component B = Component::By;

//   FDTD::FDTD field(Nx_Ny_Nz, ax_ay_az, bx_by_bz, dt);
//   Test_obj test(E, B, std::move(field));
//   Test_obj another_test(test);

//   test.set_default_field(E, B, Shift::unshifted);
//   test.analytical_default_solution(E, B, dt * t, Shift::unshifted);
//   test.numerical_solution(t, Shift::unshifted);

//   // Количество итераций 
//   uint8_t num_iterations = 6;

//   for (int64_t i{0}; i < num_iterations; ++i) {
//     dt /= 4;
//     t *= 4;
//     increase_Nx_Ny_Nz(Nx_Ny_Nz, 2); // Увеличение сетки в два раза
//     another_test = {E, B, {Nx_Ny_Nz, ax_ay_az, bx_by_bz, dt}};
//     another_test.set_default_field(E, B, Shift::unshifted);
//     another_test.analytical_default_solution(E, B, dt * t, Shift::unshifted);
//     another_test.numerical_solution(t, Shift::unshifted);
//     test.print_convergence(another_test);
//   }
// }

// TEST(Test, main_test) {
//  std::pair<int64_t, int64_t> Nx_Ny = { 64ull, 64ull };
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

//  for (int64_t y = 0ull; y < field.get_Ny(); ++y)
//    for (int64_t x = 0ull; x < field.get_Nx(); ++x) {

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
//  std::pair<int64_t, int64_t> Nx_Ny = { 64ull, 64ull };
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

//=======================================================================================================================

// TEST(Test_version_comparison, SHIFTED_Checking_the_convergence__1_iteration)
// {
//  std::tuple<int64_t, int64_t, int64_t> Nx_Ny_Nz = { 16ull, 16ull, 1ull };
//  std::tuple<double, double, double> ax_ay_az = { 0.0, 0.0, 0.0 };
//  std::tuple<double, double, double> bx_by_bz = { 1.0, 1.0, 1.0 };
//  double start = omp_get_wtime();
//
//  double dt = 1e-15;
//  // double t = 2e-14;
//  //double dx = (bx_by.first - ax_ay.first) / Nx_Ny.first;
//  double dx = (std::get<0>(bx_by_bz) - std::get<0>(ax_ay_az)) /
//  std::get<0>(Nx_Ny_Nz); dt = 0.25 * dx / C;
//  //double t = 1e-12;
//
//
//  //t = (bx_by.first - ax_ay.first) / C * 0.25;
//  int64_t t = 10; // Задание количества итераций
//
//  FDTD::FDTD field(Nx_Ny_Nz, ax_ay_az, bx_by_bz, dt);
//
//  Component E = Component::Ez;
//  Component B = Component::Bx;
//  Shift shift = Shift::shifted;
//
//  Test_obj test(E, B, std::move(field));
//  test.analytical_default_solution(E, B, t * dt, shift);  // t * dt OR t
//  test.set_default_field(E, B, shift);
//  test.numerical_solution(t, shift);
//
//  //=====Second field=====
//
//  Nx_Ny_Nz = std::make_tuple(32ull, 32ull, 1ull);
//  dt = dt / 2;
//  t *= 2;
//
//  //dx = (bx_by.first - ax_ay.first) / Nx_Ny.first;
//  //dt = 0.25 * dx / C;
//
//  FDTD::FDTD field_2(Nx_Ny_Nz, ax_ay_az, bx_by_bz, dt);
//
//  Test_obj test_2(E, B, std::move(field_2));
//  test_2.analytical_default_solution(E, B, t * dt, shift); // t * dt OR t
//  test_2.set_default_field(E, B, shift);
//  test_2.numerical_solution(t, shift);
//
//  double end = omp_get_wtime();
//  std::cout << "Time = " << end - start << std::endl;
//
//  test.print_convergence(test_2);
//}

//=======================================================================================================================

// TEST(Test_version_comparison,
// SHIFTED_Checking_the_convergence__several_iterations)
//{
//  std::cout << "======================Start several
//  operations======================\n";
//
//  std::tuple<int64_t, int64_t, int64_t> Nx_Ny_Nz = { 16ull, 16ull, 16ull }; //
//  starting grid std::tuple<double, double, double> ax_ay_az = { 0.0, 0.0, 0.0
//  }; std::tuple<double, double, double> bx_by_bz = { 1.0, 1.0, 1.0 };
//
//  auto multiply_by_two = [](const auto& element) { return element * 2; };
//
//  //double dx = (bx_by.first - ax_ay.first) / Nx_Ny.first;
//  double dx = (std::get<0>(bx_by_bz) - std::get<0>(ax_ay_az)) /
//  std::get<0>(Nx_Ny_Nz); double dt = 0.25 * dx / C; int64_t t = 10;
//
//  FDTD::FDTD field_1(Nx_Ny_Nz, ax_ay_az, bx_by_bz, dt);
//
//  Component E = Component::Ex;
//  Component B = Component::By;
//  Shift shift = Shift::shifted;
//
//  Test_obj test(E, B, std::move(field_1));
//
//  test.analytical_default_solution(E, B, t * dt, shift);
//  test.set_default_field(E, B, shift);
//  test.numerical_solution(t, shift);
//
//
//  //Nx_Ny = { Nx_Ny.first * 2, Nx_Ny.second * 2 };
//  //t *= 4;
//  //dt /= 4;
//  FDTD::FDTD field_2(Nx_Ny_Nz, ax_ay_az, bx_by_bz, dt);
//  Test_obj another_test(E, B, std::move(field_2));
//  //another_test.analytical_default_solution_OY(E, B, t * dt, shift);
//  //another_test.set_default_field_OY(E, B, shift);
//  //another_test.numerical_solution(t, shift);
//
//  //test.print_convergence(another_test);
//
//  int64_t Nx, Ny, Nz;
//  std::tie(Nx, Ny, Nz) = Nx_Ny_Nz;
//  for (uint8_t i = 0; i < 4; i++)
//  {
//    Nx *= 2; Ny *= 2;
//    if (Nz != 1ull) Nz *= 2;
//
//    if (shift == Shift::shifted)
//    {
//      t *= 2;
//      dt /= 2;
//    }
//    else
//    {
//      t *= 4;
//      dt /= 4;
//    }
//    //another_test.field.set_Nx_Ny(Nx_Ny.first, Nx_Ny.second);
//    another_test.field.set_Nx_Ny_Nz(Nx, Ny, Nz);
//    another_test.field.set_dt(dt);
//    //another_test.analytical_field.set_Nx_Ny(Nx_Ny.first, Nx_Ny.second);
//    another_test.analytical_field.set_Nx_Ny_Nz(Nx, Ny, Nz);
//    another_test.analytical_field.set_dt(dt);
//    another_test.analytical_default_solution(E, B, t * dt, shift);
//    another_test.set_default_field(E, B, shift);
//    another_test.numerical_solution(t, shift);
//    test.print_convergence(another_test);
//    //test = another_test;
//  }
//}

} // namespace gtest

#endif // !__TESTS_HPP__
