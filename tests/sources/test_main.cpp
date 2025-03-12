#include <mpi.h>
#include "gtest.h"

int main(int argc, char** argv) {

//   MPI_Init(&argc, &argv);
//   int rank, size;
//   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   MPI_Comm_size(MPI_COMM_WORLD, &size);

//   int Nx, Ny, Nz;
//   Nx = Ny = Nz = 3;
//   int matrix[3][3][3];

// int sizes[3]    = { Nx,  Ny, Nz };
// int subsizes[3] = { 1,   3, 1 };
// int starts[3]   = { 0,    0,  0 }; 

//   if (rank == 0) {
//     int counter = 0;
//     for (int i = 0; i < 3; ++i) {
//       for (int j = 0; j < 3; ++j) {
//         for (int k = 0; k < 3; ++k) {
//           matrix[i][j][k] = counter++; 
//           std::cout << matrix[i][j][k] << '\t';
//         }
//         std::cout << std::endl;
//       }
//     }
//   }
//   // Create the subarray datatype
//   MPI_Datatype matrix_type;
//   MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &matrix_type);
//   // MPI_Type_vector(Nz, Ny, Nx * Ny, MPI_INT, &matrix_type);
//   MPI_Type_commit(&matrix_type);

//   if (rank == 0) {
//     MPI_Send(matrix, 1, matrix_type, 1, 0, MPI_COMM_WORLD);
//   } else {
//     int buffer[9];
//     MPI_Recv(buffer, 9, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//     for (int i : buffer) {
//       std::cout << i << ' ';
//     }
//     std::cout << std::endl;
//   }

//   MPI_Finalize();
//   return 0;
  
  std::setlocale(LC_ALL, "ru_RU.UTF-8");
  MPI_Init(&argc, &argv);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  ::testing::InitGoogleTest(&argc, argv);
  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  if (rank != 0) {
    delete listeners.Release(listeners.default_result_printer());
  }

  int status = RUN_ALL_TESTS();

  MPI_Barrier(MPI_COMM_WORLD);

  if (status != 0) {
    std::cerr << "Tests failed!" << std::endl;
  }
  else {
    std::cout << "Tests passed!" << std::endl;
  }
  MPI_Finalize();
  
  return status;
}