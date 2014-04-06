#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
  const int PING_PONG_LIMIT = 10;

  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // We are assuming at least 2 processes for this task
  
  int partner_rank = (world_rank + 1) % 2;
//  int ping_pong_count = 0;
  
  if (world_rank % 2 == 0) {
      // Increment the ping pong count before you send it
      char str[20] = "Hello";
      MPI_Send(str, 6, MPI_CHAR, partner_rank, 0, MPI_COMM_WORLD);
      printf(" sent %s \n", str);
    
  } else {
      char ctr[20];
      MPI_Recv(ctr, 6, MPI_CHAR, partner_rank, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      printf(" received %s", ctr);
    }
  
  MPI_Finalize();

}
