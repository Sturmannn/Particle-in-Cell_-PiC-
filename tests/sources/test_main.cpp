#include <mpi.h>
#include "gtest.h"

int main(int argc, char** argv) {
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