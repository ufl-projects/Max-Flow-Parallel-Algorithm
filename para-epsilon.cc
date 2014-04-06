#include <iostream>
#include <fstream>
#include "mpi.h"

using namespace std;

int main(int argc, char ** argv) {
    
    int rank, size;
    char name[80];
    int length; 
    int nodes, edges;
    int source = 0;
    double **capacity; // = new float[n]();
    double **flow;

    MPI_Init(&argc, &argv); // note that argc and argv are passed
                            // by address

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Get_processor_name(name,&length);


    if(rank==0){	// master, read graph;
   
	int n1,n2,n3;

	ifstream f;	// read network graph
	f.open("graph.txt");
	if(!f){
		 cout<<"Errror opening file";
		return 0;
	}

	f >> nodes >> edges;

	capacity = new double*[nodes+1];
	flow = new double*[nodes+1];

	for(int i = 0; i < nodes; ++i){
	    capacity[i] = new double[nodes+1]();
	    flow[i] = new double[nodes+1]();
	}

	while(!f.eof())
	{

		f >> n1 >> n2 >> n3;

		capacity[n1][n2] = n3;
	}

	cout << "file read";

	for(int i=1;i<nodes;i++)		// preflow operation
		for(int j=1;j<nodes;j++)
			if(i==source) 
				flow[i][j] = capacity[i][j];

	// broadcast to other nodes

	// send number of nodes

	// send capacity matrix

	// send flow matrix

	f.close();
    }

    else{	// slave nodes

	// receive nodes, edges

	// receive capacity matrix

	// receive flow matrix

	// generate push list

	// get edge from push list

	//  

	


    }



   
    MPI_Finalize();
}


/* Not added

void push(const int * const * C, int ** F, int *excess, int u, int v) {
  int send = MIN(excess[u], C[u][v] - F[u][v]);
  F[u][v] += send;
  F[v][u] -= send;
  excess[u] -= send;
  excess[v] += send;
}
 
void relabel(const int * const * C, const int * const * F, int *height, int u) {
  int v;
  int min_height = INFINITE;
  for (v = 0; v < NODES; v++) {
    if (C[u][v] - F[u][v] > 0) {
      min_height = MIN(min_height, height[v]);
      height[u] = min_height + 1;
    }
  }
};
 
void discharge(const int * const * C, int ** F, int *excess, int *height, int *seen, int u) {
  while (excess[u] > 0) {
    if (seen[u] < NODES) {
      int v = seen[u];
      if ((C[u][v] - F[u][v] > 0) && (height[u] > height[v])){
	push(C, F, excess, u, v);
      }
      else
	seen[u] += 1;
    } else {
      relabel(C, F, height, u);
      seen[u] = 0;
    }
  }
}
 
void moveToFront(int i, int *A) {
  int temp = A[i];
  int n;
  for (n = i; n > 0; n--){
    A[n] = A[n-1];
  }
  A[0] = temp;
}
 
int pushRelabel(const int * const * C, int ** F, int source, int sink) {
  int *excess, *height, *list, *seen, i, p;
 
  excess = (int *) calloc(NODES, sizeof(int));
  height = (int *) calloc(NODES, sizeof(int));
  seen = (int *) calloc(NODES, sizeof(int));
 
  list = (int *) calloc((NODES-2), sizeof(int));
 
  for (i = 0, p = 0; i < NODES; i++){
    if((i != source) && (i != sink)) {
      list[p] = i;
      p++;
    }
  }
 
  height[source] = NODES;
  excess[source] = INFINITE;
  for (i = 0; i < NODES; i++)
    push(C, F, excess, source, i);
 
  p = 0;
  while (p < NODES - 2) {
    int u = list[p];
    int old_height = height[u];
    discharge(C, F, excess, height, seen, u);
    if (height[u] > old_height) {
      moveToFront(p,list);
      p=0;
    }
    else
      p += 1;
  }
  int maxflow = 0;
  for (i = 0; i < NODES; i++)
    maxflow += F[source][i];
 
  free(list);
 
  free(seen);
  free(height);
  free(excess);
 
  return maxflow;
}

*/
