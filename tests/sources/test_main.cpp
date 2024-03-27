#include <mpi.h>
#include "gtest.h"

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  ::testing::InitGoogleTest(&argc, argv);
  int status = RUN_ALL_TESTS();
  MPI_Finalize();
  return status;
}