#include <mpi.h>
#include "gtest.h"
#include <set>

int main(int argc, char** argv) {

  std::set<int> set;
  set.insert(2)
  // const int matrix_size = 4;

  // int matrix[matrix_size][matrix_size][matrix_size];
  // MPI_Init(&argc, &argv);
  // int rank, size;
  // MPI_Comm_size(MPI_COMM_WORLD, &size);
  // MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // if (rank == 0) {
  //   int counter = 0;
  //   for (int i = 0; i < matrix_size; ++i) {
  //     for (int j = 0; j < matrix_size; ++j) {
  //       for (int k = 0; k < matrix_size; ++k) {
  //         matrix[i][j][k] = counter++;
  //         std::cout << matrix[i][j][k] << '\t'; 
  //       }
  //       std::cout << std::endl;
  //     }
  //     std::cout << std::endl;
  //   }
  // }

  // MPI_Datatype col;
  // MPI_Type_vector(matrix_size, matrix_size, matrix_size * matrix_size, MPI_INT, &col);
  // MPI_Type_commit(&col);
  
  // int received_face[matrix_size][matrix_size];
  // if (rank == 0) {
  //   MPI_Send(&(matrix[0][0][0]), 1, col, 1, 0, MPI_COMM_WORLD);
  // } 
  // else if (rank == 1) {
  //   MPI_Recv(&(received_face[0][0]), matrix_size * matrix_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //   for (int i = 0; i < matrix_size; ++i) {
  //     for (int j = 0; j < matrix_size; ++j)  {
  //       std::cout << received_face[i][j] << ' ';
  //     }
  //     std::cout << std::endl;
  //   }
  //   std::cout << std::endl;
  // }

  // MPI_Finalize();
  // return 0;


  
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