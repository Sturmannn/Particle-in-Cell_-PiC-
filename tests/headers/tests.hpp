#ifndef __TESTS_HPP__
#define __TESTS_HPP__

#include <functional>
#include "FDTD.hpp"
#include "gtest.h"


//using FDTD::Axis;
using FDTD::Component;

namespace gtest {

enum class Shift { shifted, unshifted };

class Test_obj {
public:
  Test_obj() = delete;
  Test_obj(const Component _E, const Component _B, int64_t Nx, int64_t Ny); // Вариант пока что для 2D
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
  
  int (&get_coords(void))[2] { return coords; }
  int (&get_dims(void))[2] { return dims; }

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
  int dims[2], coords[2]; // Размеры топологии и координаты процесса
  int64_t local_Nx, local_Ny, local_Nz; // Размеры поддомена

  void Courant_condition_check(const Shift _shift) const noexcept;
  void set_default_field_or_analytical_default_solution(const Component E,
                                                        const Component B,
                                                        const Shift _shift,
                                                        std::function<void(std::tuple<Axis, int64_t, int64_t>,
                                                          std::tuple<Axis, int64_t, int64_t>,
                                                          std::tuple<Axis, int64_t, int64_t>)> loop_function);
  void create_cartesian_topology(int world_size);
  
  // Пока что здесь только 2D случай
  void set_subdomain_sizes(int64_t Nx, int64_t Ny);
};

// Пример декартовой топологии для 2D сетки (но без периодических граничных условий)
// [ 6 ]---[ 7 ]---[ 8 ]
//   |       |       |
// [ 3 ]---[ 4 ]---[ 5 ]
//   |       |       |
// [ 0 ]---[ 1 ]---[ 2 ]

// Декартова топология для MPI
// void create_cartesian_topology(int world_size, int dims[2], int coords[2], MPI_Comm &cart_comm) {
//     // Определение размеров декартовой топологии (сколько процессов вдоль каждой из осей)
//     dims[0] = dims[1] = 0;
//     MPI_Dims_create(world_size, 2, dims); // 2 - количество измерений (для 2D сетки)

//     // Создание декартовой топологии
//     int periods[2] = {1, 1}; // Для периодических граничных условий
//     // Возможно с переодичностью будет попроще, но не факт.

//     MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);

//     // Получение координаты текущего процесса в топологии
//     int rank;
//     MPI_Comm_rank(cart_comm, &rank);
//     MPI_Cart_coords(cart_comm, rank, 2, coords);
// }

// void set_subdomain_sizes(int64_t Nx, int64_t Ny, int dims[2], int coords[2], int64_t &local_Nx, int64_t &local_Ny) {
//     // Пока что здесь только 2D случай

//     // Определение размера блоков матрицы для каждого процесса
//     local_Nx = (Nx + 2 * dims[0]) / dims[0]; // +2 для граничных полей
//     local_Ny = (Ny + 2 * dims[1]) / dims[1]; // +2 для граничных полей

//     // Учет остатка, если количество размер сетки не кратен количеству процессов
//     if (coords[0] < (Nx + 2 * dims[0]) % dims[0]) local_Nx++;
//     if (coords[1] < (Ny + 2 * dims[1]) % dims[1]) local_Ny++;

//     // Уменьшение на 2, так как в конструкторе учитываются граничные поля
//     local_Nx -= 2;
//     local_Ny -= 2;
// }

TEST(Test_version_comparison, shifted_OZ) {
  // Потенциальны проблемы при запуске:
  // 1. Ось по OZ, а сетка 2D
  std::tuple<int64_t, int64_t, int64_t> Nx_Ny_Nz = {10, 10, 1};
  std::tuple<double, double, double> ax_ay_az = {0.0, 0.0, 0.0};
  std::tuple<double, double, double> bx_by_bz = {1.0, 1.0, 1.0};
  std::tuple<double, double, double> dx_dy_dz = 
        {(std::get<0>(bx_by_bz) - std::get<0>(ax_ay_az)) / static_cast<double>(std::get<0>(Nx_Ny_Nz)),
         (std::get<1>(bx_by_bz) - std::get<1>(ax_ay_az)) / static_cast<double>(std::get<1>(Nx_Ny_Nz)),
         (std::get<2>(bx_by_bz) - std::get<2>(ax_ay_az)) / static_cast<double>(std::get<2>(Nx_Ny_Nz))};
  
  if (std::get<2>(Nx_Ny_Nz) > 1) {
    std::get<2>(dx_dy_dz) = 0.0;
  }
  // double dt = 2e-15;

  double dx = std::get<0>(dx_dy_dz);
  double dt = 0.25 * dx / C;
  int64_t t = 250; // Задание количества итераций

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // int dims[2], coords[2];
    // MPI_Comm cart_comm;
    // create_cartesian_topology(world_size, dims, coords, cart_comm);

  Component E = Component::Ey;
  Component B = Component::Bz;
  Shift shift = Shift::shifted;
  
  if (rank == 0) std::cout << "Axis: " << FDTD::FDTD::axisToString(E, B) << std::endl;



// =========================================================

// Пока что MPI для 2D и разбиение по процессам идёт по оси OY
// В дальнейшем скорее всего должен быть выбор, по какой из осей разбивать
// НЕ ЗАБЫТЬ ТО ЖЕ СДЕЛАТЬ ДЛЯ ВТОРОГО ТЕСТОВОГО ПРИМЕРА

  // int64_t local_Nx = std::get<0>(Nx_Ny_Nz);
  // int64_t local_Ny = (std::get<1>(Nx_Ny_Nz) + 2 * world_size) / world_size;

  // // Это рассчитывается с учётом дополнительных граничных полей.
  // // В реализации же этого кода, они учитываются в конструктуре, поэтому этот результат нужно будет уменьшить на 2
  // if (rank < (std::get<1>(Nx_Ny_Nz) + 2 * world_size) % world_size) local_Ny += 1;
  
  // // Уменьшаем на 2, так как в конструкторе учитываются граничные поля
  // // if (world_size > 1) local_Ny -= 2;
  // local_Ny -= 2;


  // int64_t local_Nx, local_Ny;
  // set_subdomain_sizes(std::get<0>(Nx_Ny_Nz), std::get<1>(Nx_Ny_Nz), dims, coords, local_Nx, local_Ny);

  // std::cout << "Rank = " << rank << " Local_Nx = " << local_Nx << " Local_Ny = " << local_Ny << std::endl;

  // Пока что Nz = 1
  // std::tuple<int64_t, int64_t, int64_t> local_Nx_Ny_Nz = {local_Nx, local_Ny, 1};


// =========================================================

  Test_obj test(E, B, std::get<0>(Nx_Ny_Nz), std::get<1>(Nx_Ny_Nz)); // Пока что 2D
  // Внимательно проверять, что в field нужно передавать именно local размеры
  FDTD::FDTD field(test.get_local_Nx_Ny_Nz(), ax_ay_az, bx_by_bz, dx_dy_dz, dt); 
  test.initialize_field(field);


  // FDTD::FDTD field(local_Nx_Ny_Nz, ax_ay_az, bx_by_bz, dx_dy_dz, dt);
  // Test_obj test(E, B, std::move(field));

  // std::vector<double>& vector = test.get_field().get_Ex().get_field();
  // if (rank == 0) {
  //   for (int64_t i = 0; i < vector.size(); ++i) {
  //     vector[i] = i;
  //     std::cout << vector[i] << " ";
  //   }
  //   std::cout << '\n' << std::endl;
  // }
  //
  // MPI_Barrier(MPI_COMM_WORLD);
  // if (rank == 1) {
  //   for (int64_t i = 0; i < vector.size(); ++i) {
  //     vector[i] = i + 48;
  //     std::cout << vector[i] << " ";
  //   }
  //   std::cout << '\n' << std::endl;
  // }
  // MPI_Barrier(MPI_COMM_WORLD);
  // if (rank == 2) {
  //   for (int64_t i = 0; i < vector.size(); ++i) {
  //     vector[i] = i + 84;
  //     std::cout << vector[i] << " ";
  //   }
  //   std::cout << '\n' << std::endl;
  // }
  // MPI_Barrier(MPI_COMM_WORLD);
  //
  //
  // test.get_field().boundary_synchronization();
  // if (rank == 0) {
  //   for (int64_t i = 0; i < vector.size(); ++i) {
  //     std::cout << test.get_field().get_Ex().get_field()[i] << " ";
  //   }
  //   std::cout << '\n' << std::endl;
  // }
  //
  // MPI_Barrier(MPI_COMM_WORLD);
  // if (rank == 1) {
  //   for (int64_t i = 0; i < vector.size(); ++i) {
  //     std::cout << test.get_field().get_Ex().get_field()[i] << " ";
  //   }
  //   std::cout << '\n' << std::endl;
  // }
  // MPI_Barrier(MPI_COMM_WORLD);
  // if (rank == 2) {
  //   for (int64_t i = 0; i < vector.size(); ++i) {
  //     std::cout << test.get_field().get_Ex().get_field()[i] << " ";
  //   }
  //   std::cout << '\n' << std::endl;
  // }
  // MPI_Barrier(MPI_COMM_WORLD);
  //
  // exit(-1);

  test.analytical_default_solution(E, B, t * dt, shift);
  test.set_default_field(E, B, shift);
  test.numerical_solution(t, shift, test.get_cart_comm());


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
