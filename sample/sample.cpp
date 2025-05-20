#include <chrono>
#include "FDTD_MPI.hpp"
#include "FieldFileManager.hpp"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank, world_size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  const int Nx = 64;
  const int Ny = 64;
  const int Nz = 64;

  FDTD::UnboundedGridSizes grid_sizes = {Nx, Ny, Nz};
  FDTD::GridCoordinatesBounds bounds = {1.0, 2.0, 3.0, 2.0,
                                        3.0, 4.0}; // ax, ay, az, bx, by, bz

  const double dt = 0.25 * (bounds.bx - bounds.ax) /
              static_cast<double>(grid_sizes.Nx) / FDTD::C;
  const int iterations = 150;
  
  if (rank == 0) {
    std::cout << "OMP max threads: " << omp_get_max_threads() << std::endl;
  }
  // #pragma omp parallel
  // {
  //   if (rank == 0) {
  //     std::cout << "Number of thread: " << omp_get_thread_num() << std::endl;
  //   }
  // }
  FDTD::Component E = FDTD::Component::Ey;
  FDTD::Component B = FDTD::Component::Bz;
  // FDTD::Component E = FDTD::Component::Ex;
  // FDTD::Component B = FDTD::Component::By;
  // FDTD::Component E = FDTD::Component::Ex;
  // FDTD::Component B = FDTD::Component::Bz;
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
    // Writing time to file -------
    std::array<int, 3> dims = mpi_wrapper_ptr->get_dims();
    std::ofstream measure_file;
    measure_file.open(FDTD::path_to_measurements_file, std::ios::app);
    if (measure_file.is_open()) {
      measure_file << elapsed.count() << '\t' << omp_get_max_threads() << '\t' << world_size << "\t" <<
        dims[0] << '\t' << dims[1] << '\t' << dims[2] << '\t' << Nx  << '\t' <<
          Ny << '\t' << Nz << std::endl;
    }
    else {
      std::cerr << "The file for measurements can't be opened" << std::endl;
    }
    measure_file.close();
    
    // ----------------------------
  }


  FDTD::FieldFileManager field_file_manager(mpi_wrapper_ptr);
  field_file_manager.clear_files(FDTD::path_to_analytical_data_directory);
  field_file_manager.write_fields_to_file(
      FDTD::path_to_analytical_data_directory, E, B,
      analytical_solver.get_space_delta(E, B), analytical_solver, 0);

  field_file_manager.clear_files(FDTD::path_to_calculated_data_directory);
  field_file_manager.write_fields_to_file(
      FDTD::path_to_calculated_data_directory, E, B,
      numerical_solver.get_space_delta(E, B), numerical_solver, 0);
  mpi_wrapper_ptr->finalize();
  MPI_Finalize();
  return 0;
}