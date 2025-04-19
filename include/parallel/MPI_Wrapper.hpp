#include "../сore/Grid.hpp"
#include "../сore/types.hpp"
#include <array>
#include <memory>
#include <mpi.h>

namespace FDTD {

class MPI_Wrapper {
public:
  MPI_Wrapper(MPI_Comm _world_comm, Grid &_grid);
  ~MPI_Wrapper() {
    if (cart_comm != MPI_COMM_NULL)
      MPI_Comm_free(&cart_comm);
  };

  int get_world_size() const { return world_size; }
  std::array<int, 3> get_coords() const { return coords; }
  int get_cart_rank(int x, int y, int z) const {
    int rank = 0;
    std::array<int, 3> coords = {x, y, z};
    MPI_Cart_rank(cart_comm, coords.data(), &rank);
    return rank;
  }

  std::vector<int> get_sizes_OX() const {
    std::vector<int> all_local_Nx(dims[0], 0);
  
    MPI_Allgather(grid->get_mpi_local_unbounded_subdomain_sizes().Nx, 1, MPI_INT,
                  all_local_Nx.data(), 1, MPI_INT, world_comm);
    return all_local_Nx;
  }
  std::vector<int> get_sizes_OY() const {
    std::vector<int> all_local_Ny(dims[1], 0);
    MPI_Allgather(grid->get_mpi_local_unbounded_subdomain_sizes().Ny, 1, MPI_INT,
                  all_local_Ny.data(), 1, MPI_INT, world_comm);
    return all_local_Ny;
  }
  std::vector<int> get_sizes_OZ() const {
    std::vector<int> all_local_Nz(dims[2], 0);
    MPI_Allgather(grid->get_mpi_local_unbounded_subdomain_sizes().Nz, 1, MPI_INT,
                  all_local_Nz.data(), 1, MPI_INT, world_comm);
    return all_local_Nz;
  }

private:
  void create_cartesian_topology(int world_size);
  void set_subdomain_sizes(void);

  std::unique_ptr<Grid> grid;
  std::array<int, 3> dims;
  std::array<int, 3> coords;
  MPI_Comm world_comm;
  MPI_Comm cart_comm;
  int world_size;
  int world_rank;
};

inline FDTD::MPI_Wrapper::MPI_Wrapper(MPI_Comm _world_comm, Grid &_grid)
    : grid(std::make_unique<Grid>(_grid)), world_comm(_world_comm) {
  MPI_Comm_size(world_comm, &world_size);
  MPI_Comm_rank(world_comm, &world_rank);

  create_cartesian_topology(world_size);

  UnboundedGridSizes unbounded_grid_sizes = grid->get_unbounded_grid_sizes();

  // int mpi_local_Nx = (unbounded_grid_sizes.Nx + )
}
inline void MPI_Wrapper::create_cartesian_topology(int world_size) {
  // Periodic boundary conditions
  std::array<int, MPI_DIMENSION> periods = {1, 1, 1};

  // Define the dimensions of the Cartesian topology
  MPI_Dims_create(world_size, MPI_DIMENSION, dims.data());
  // Create the Cartesian topology
  MPI_Cart_create(world_comm, MPI_DIMENSION, dims.data(), periods.data(), 1,
                  &cart_comm);
  // Get the rank of the current process in the Cartesian topology
  MPI_Cart_coords(cart_comm, world_rank, MPI_DIMENSION, coords.data());

  if (world_rank == 0) {
    std::cout << "MPI dimension: " << MPI_DIMENSION << "D" << std::endl;
    std::cout << "Count of processes: " << world_size << std::endl;
    std::cout << "Ox processes: " << dims[0] << " Oy processes: " << dims[1]
              << " Oz processes: " << dims[2] << std::endl;
    std::cout << std::endl;
  }
}

inline void MPI_Wrapper::set_subdomain_sizes(void) {
  UnboundedGridSizes unbounded_grid_sizes = grid->get_unbounded_grid_sizes();
  BoundedGridSizes mpi_local_bounded_subdomain_sizes;


  mpi_local_bounded_subdomain_sizes.Nx =
      (unbounded_grid_sizes.Nx + 2 * dims[0]) /
      dims[0]; // +2 for boundary fields
  mpi_local_bounded_subdomain_sizes.Ny =
      (unbounded_grid_sizes.Ny + 2 * dims[1]) /
      dims[1]; // +2 for boundary fields
  mpi_local_bounded_subdomain_sizes.Nz =
      (unbounded_grid_sizes.Nz + 2 * dims[2]) /
      dims[2]; // +2 for boundary fields

  if (coords[0] < (unbounded_grid_sizes.Nx + 2 * dims[0]) % dims[0])
    mpi_local_bounded_subdomain_sizes.Nx++;
  if (coords[1] < (unbounded_grid_sizes.Ny + 2 * dims[1]) % dims[1])
    mpi_local_bounded_subdomain_sizes.Ny++;
  if (coords[2] < (unbounded_grid_sizes.Nz + 2 * dims[2]) % dims[2])
    mpi_local_bounded_subdomain_sizes.Nz++;

  grid->set_mpi_local_bounded_subdomain_sizes(
      mpi_local_bounded_subdomain_sizes);

  grid->set_mpi_local_unbounded_subdomain_sizes({
    mpi_local_bounded_subdomain_sizes.Nx - 2,
    mpi_local_bounded_subdomain_sizes.Ny - 2,
    mpi_local_bounded_subdomain_sizes.Nz - 2
  });
}

} // namespace FDTD