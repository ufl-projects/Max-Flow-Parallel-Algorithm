#include <iostream>
#include "mpi.h"

using namespace std;

int main(int argc, char ** argv) {
    int rank, size;
    char name[80];
    int length; 

    MPI_Init(&argc, &argv); // note that argc and argv are passed
                            // by address

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Get_processor_name(name,&length);

    
 
   int number;
   if (rank == 0) {
     number = -1;
     MPI_Send(&number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
   } else if (rank == 1) {
    MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    cout << "Process 1 received number " << number << " from process 0\n";
    
   }



    cout<< "Hello MPI: processor " << rank << " of " << size << "on " << name << "\n";
   
    MPI_Finalize();
}
