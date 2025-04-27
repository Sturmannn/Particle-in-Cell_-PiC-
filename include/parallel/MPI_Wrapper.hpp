#include "Grid.hpp"
#include "types.hpp"
#include <array>
#include <memory>
#include <mpi.h>

namespace FDTD {

class MPI_Wrapper {
public:
  MPI_Wrapper(MPI_Comm _world_comm, std::shared_ptr<Grid> _grid);
  ~MPI_Wrapper() {
    if (cart_comm != MPI_COMM_NULL)
      MPI_Comm_free(&cart_comm);
  };

  int get_world_size() const { return world_size; }
  int get_world_rank() const { return world_rank; }
  MPI_Comm get_cart_comm() const { return cart_comm; }
  std::array<int, 3> get_coords() const { return coords; }

  int get_cart_rank(int x, int y, int z) const {
    int rank = 0;
    std::array<int, 3> coords = {x, y, z};
    MPI_Cart_rank(cart_comm, coords.data(), &rank);
    return rank;
  }

  std::vector<int> get_sizes_OX() const {
    std::vector<int> all_local_Nx(world_size, 0);
    int Nx = grid->get_mpi_local_unbounded_subdomain_sizes().Nx;
    MPI_Allgather(&Nx, 1, MPI_INT, all_local_Nx.data(), 1, MPI_INT, world_comm);
    return all_local_Nx;
  }
  std::vector<int> get_sizes_OY() const {
    std::vector<int> all_local_Ny(world_size, 0);
    int Ny = grid->get_mpi_local_unbounded_subdomain_sizes().Ny;
    MPI_Allgather(&Ny, 1, MPI_INT, all_local_Ny.data(), 1, MPI_INT, world_comm);
    return all_local_Ny;
  }
  std::vector<int> get_sizes_OZ() const {
    std::vector<int> all_local_Nz(world_size, 0);
    int Nz = grid->get_mpi_local_unbounded_subdomain_sizes().Nz;
    MPI_Allgather(&Nz, 1, MPI_INT, all_local_Nz.data(), 1, MPI_INT, world_comm);
    return all_local_Nz;
  }

  void mpi_sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag,
                    void *recvbuf, int recvcount, MPI_Datatype recvtype, int source,
                    int recvtag, MPI_Comm comm, MPI_Status *status) {
    MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount,
                 recvtype, source, recvtag, comm, status);
  }

private:
  void create_cartesian_topology(int world_size);
  void set_subdomain_sizes(void);

  std::shared_ptr<Grid> grid;
  std::array<int, 3> dims;
  std::array<int, 3> coords;
  MPI_Comm world_comm;
  MPI_Comm cart_comm;
  int world_size;
  int world_rank;
};

inline FDTD::MPI_Wrapper::MPI_Wrapper(MPI_Comm _world_comm,
                                      std::shared_ptr<Grid> _grid)
    : grid(_grid), world_comm(_world_comm) {
  MPI_Comm_size(world_comm, &world_size);
  MPI_Comm_rank(world_comm, &world_rank);

  dims.fill(0);
  coords.fill(0);

  create_cartesian_topology(world_size);

  UnboundedGridSizes unbounded_grid_sizes = grid->get_unbounded_grid_sizes();

  set_subdomain_sizes();
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

  grid->set_mpi_local_unbounded_subdomain_sizes(
      {mpi_local_bounded_subdomain_sizes.Nx - 2,
       mpi_local_bounded_subdomain_sizes.Ny - 2,
       mpi_local_bounded_subdomain_sizes.Nz - 2});
}

} // namespace FDTD