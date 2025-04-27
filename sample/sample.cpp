#include "FDTD_MPI.hpp"
#include "FieldFileManager.hpp"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    FDTD::UnboundedGridSizes grid_sizes = {16, 16, 16}; // Nx, Ny, Nz
    FDTD::GridCoordinatesBounds bounds = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0}; // ax, ay, az, bx, by, bz
    
    double dt = 0.25 * (bounds.bx - bounds.ax) / static_cast<double>(grid_sizes.Nx) / FDTD::C;
    int iterations = 150;

    FDTD::Component E = FDTD::Component::Ez;
    FDTD::Component B = FDTD::Component::Bx;
    FDTD::Shift shift = FDTD::Shift::shifted;

    auto grid_ptr = std::make_shared<FDTD::Grid>(grid_sizes, bounds, dt);
    auto mpi_wrapper_ptr = std::make_shared<FDTD::MPI_Wrapper>(MPI_COMM_WORLD, grid_ptr);


    FDTD::AnalyticalSolverFDTD analytical_solver(grid_ptr, mpi_wrapper_ptr);
    analytical_solver.solve(E, B, iterations * dt, shift); // t = 0
    std::cout << "Rank: " << rank << " Hellooo!" << std::endl;

    FDTD::NumericalSolverFDTD numerical_solver(grid_ptr, mpi_wrapper_ptr);
    numerical_solver.solve(E, B, iterations, shift);
    // FDTD::GridCoordinatesSteps steps = grid_ptr->get_coordinates_steps(); // dx, dy, dz
    

    FDTD::FieldFileManager field_file_manager(mpi_wrapper_ptr);
    field_file_manager.clear_files(FDTD::path_to_analytical_data_directory);
    field_file_manager.write_fields_to_file(FDTD::path_to_analytical_data_directory, E, B, analytical_solver.get_space_delta(E, B), analytical_solver, 0);
    
    field_file_manager.clear_files(FDTD::path_to_calculated_data_directory);
    field_file_manager.write_fields_to_file(FDTD::path_to_calculated_data_directory, E, B, numerical_solver.get_space_delta(E, B), numerical_solver, 0);

    mpi_wrapper_ptr->~MPI_Wrapper();
    MPI_Finalize();
    return 0;
}