#include <mpi.h>
#include "gtest.h"

int main(int argc, char** argv) {
  std::setlocale(LC_ALL, "ru_RU.UTF-8");
  MPI_Init(&argc, &argv);

  ::testing::InitGoogleTest(&argc, argv);
  int status = RUN_ALL_TESTS();
  MPI_Finalize();
  
  if (status != 0) {
    std::cerr << "Tests failed!" << std::endl;
  }
  else {
    std::cout << "Tests passed!" << std::endl;
  }
  
  return status;

}