#include <mpi.h>
#include "gtest.h"
// #include "tests.hpp"
#include <set>

namespace gtest {
  // Default values
  int arg_Nx = 16;
  int arg_Ny = 16;
  int arg_Nz = 16;
}

int main(int argc, char** argv) {


  // const int matrix_size = 4;

  // int matrix[matrix_size][matrix_size][matrix_size];
  // std::vector<int> buffer(matrix_size * matrix_size * matrix_size, -1);

  // MPI_Init(&argc, &argv);
  // int rank, size;
  // MPI_Comm_size(MPI_COMM_WORLD, &size);
  // MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // if (rank == 0) {
  //   int counter = 0;
  //   for (int i = 0; i < matrix_size * matrix_size * matrix_size; ++i) {
  //     buffer[i] = i;
  //   }
    
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

  // MPI_Datatype subarray_send;
  // int array_of_sizes[] = {4, 4, 4};
  // int array_of_subsizes[] = {4, 2, 1};
  // int array_of_start[] = {0, 1, 3};

  // int received_face[matrix_size][matrix_size];

  // MPI_Type_create_subarray(3, array_of_sizes, array_of_subsizes, array_of_start,
  // MPI_ORDER_C, MPI_INT, &subarray_send);
  // MPI_Type_commit(&subarray_send);

  // // int recv_subarray[4][4];
  // MPI_Datatype subarray_recv;
  // MPI_Type_create_subarray(3, array_of_sizes, array_of_subsizes, array_of_start,
  //   MPI_ORDER_C, MPI_INT, &subarray_recv);

  // MPI_Type_commit(&subarray_recv);
  // int recv_matrix[matrix_size][matrix_size][matrix_size];
  // if (rank == 0) {
  //   MPI_Send(&matrix[0][0][0], 1, subarray_send, 1, 0, MPI_COMM_WORLD);
  // }
  // else if (rank == 1) {
  //   // MPI_Recv(&(recv_subarray[0][0]), matrix_size * matrix_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  //   int count = 0;
  //   MPI_Recv(buffer.data(), 1, subarray_send, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //   for (int i = 0; i < matrix_size * matrix_size * matrix_size; ++i) {
  //     if (i % 4 == 0 && i != 0) std::cout << "\n";
  //     std::cout << buffer[i] << ' ';
  //   }

  //   for (int i = 0; i < matrix_size; ++i) {
  //     for (int j = 0; j < matrix_size; ++j)  {
  //       for (int k = 0; k < matrix_size; ++k) {
  //         std::cout << recv_matrix[i][j][k] << ' ';
  //       }
  //       std::cout << std::endl;
  //     }
  //     std::cout << std::endl;
  //   }
  // }
  // MPI_Finalize();
  // return 0;

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
  
  // Обработка аргументов командной строки
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--Nx" && i + 1 < argc) {
      gtest::arg_Nx = std::atoi(argv[++i]);
    } else if (arg == "--Ny" && i + 1 < argc) {
      gtest::arg_Ny = std::atoi(argv[++i]);
    } else if (arg == "--Nz" && i + 1 < argc) {
      gtest::arg_Nz = std::atoi(argv[++i]);
    } else {
      std::cerr << "Command's argument ERROR\n";
    }
    
  }

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