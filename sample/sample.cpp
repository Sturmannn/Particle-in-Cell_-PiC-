#include <chrono>
#include "FDTD_MPI.hpp"
#include "FieldFileManager.hpp"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  FDTD::UnboundedGridSizes grid_sizes = {64, 64, 64}; // Nx, Ny, Nz
  FDTD::GridCoordinatesBounds bounds = {0.0, 0.0, 0.0, 1.0,
                                        1.0, 1.0}; // ax, ay, az, bx, by, bz

  double dt = 0.25 * (bounds.bx - bounds.ax) /
              static_cast<double>(grid_sizes.Nx) / FDTD::C;
  int iterations = 150;

  if (rank == 0) {
    std::cout << "OMP max threads: " << omp_get_max_threads() << std::endl;
  }

  FDTD::Component E = FDTD::Component::Ey;
  FDTD::Component B = FDTD::Component::Bz;
  FDTD::Shift shift = FDTD::Shift::shifted;

  auto grid_ptr = std::make_shared<FDTD::Grid>(grid_sizes, bounds, dt);
  auto mpi_wrapper_ptr =
      std::make_shared<FDTD::MPI_Wrapper>(MPI_COMM_WORLD, grid_ptr);

  FDTD::AnalyticalSolverFDTD analytical_solver(grid_ptr, mpi_wrapper_ptr);
  analytical_solver.solve(E, B, iterations * dt, shift); // t = 0
  
  FDTD::NumericalSolverFDTD numerical_solver(grid_ptr, mpi_wrapper_ptr);
  auto start = std::chrono::high_resolution_clock::now();
  numerical_solver.solve(E, B, iterations, shift);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  if (rank == 0) {
    std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;
  }

  FDTD::FieldFileManager field_file_manager(mpi_wrapper_ptr);
  field_file_manager.clear_files(FDTD::path_to_analytical_data_directory);
  field_file_manager.write_fields_to_file(
      FDTD::path_to_analytical_data_directory, E, B,
      analytical_solver.get_space_delta(E, B), analytical_solver, 0);
  std::cout << "Analytical solution written to file." << std::endl;

  field_file_manager.clear_files(FDTD::path_to_calculated_data_directory);
  field_file_manager.write_fields_to_file(
      FDTD::path_to_calculated_data_directory, E, B,
      numerical_solver.get_space_delta(E, B), numerical_solver, 0);
  std::cout << "Numerical solution written to file." << std::endl;
  mpi_wrapper_ptr->finalize();
  mpi_wrapper_ptr.reset();
  grid_ptr.reset();
  std::cout << "MPI_Wrapper destructor called." << std::endl;
  MPI_Finalize();
  std::cout << "MPI_Finalize called." << std::endl;
  return 0;
}