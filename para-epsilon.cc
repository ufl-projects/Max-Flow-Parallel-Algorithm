#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include "mpi.h"


#define NODES 6
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define INFINITE 10000000

using namespace std;

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

void printMatrix(const double * const * M) {
    int i,j;
    for (i = 0; i < NODES; i++) {
        for (j = 0; j < NODES; j++)
            printf("%f\t",M[i][j]);
        printf("\n");
    }
}

int main(void) {

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    string capacitystring;
    string flowstring;

    double **flow, **capacity,*excess;
    int i,nodes,edges;
    int source = 0;
    int sink = 5;
    nodes=NODES;
	


    capacity = (double **) calloc(NODES, sizeof(double*));
    for (i = 0; i < NODES; i++) {
        flow[i] = (double *) calloc(NODES, sizeof(double));
        capacity[i] = (double *) calloc(NODES, sizeof(double));
    }



    //Slave process
    if(rank!=0){

        double excess=0,push_flow=0;

        MPI_Bcast(&nodes,1, MPI_INT, rank, MPI_COMM_WORLD);
	cout <<  "My nodes " << nodes << "\n";
        MPI_Bcast(&(capacity[0][0]),nodes*nodes, MPI_DOUBLE, rank, MPI_COMM_WORLD);
        MPI_Bcast(&(flow[0][0]), nodes*nodes, MPI_DOUBLE, rank, MPI_COMM_WORLD);

	cout << rank << ": Received Capacity and Flow matrices \n";        

	//cout << "Capacity \n";
	//printMatrix(capacity);
	//printMatrix(flow);
	double* slaveFlow = new double[nodes+1]();
	
        //Receive excess from adjacent nodes
	//cout << "Outside for";
        for(int i = 0; i<rank; i++)
        {
	    if( capacity[i][rank] != 0 ) {
                double t;
		MPI_Recv(&t,1,MPI_DOUBLE, i, 2, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                slaveFlow[i] = t;
                excess += t;
		cout << rank << ": Received Initial Flow from " << i << " = " << t << "\n";     
            }
        }
	
        /*for(int i = rank+1; i<sink; i++)
        {
            if ( capacity[rank][i] != 0 && excess != 0 )
            {
                push_flow = MIN( capacity[rank][i] , excess );
                slaveFlow[i] = push_flow;
                excess -= push_flow;
                MPI_Send(&push_flow, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
            }
        }*/
	/*

        double sendExcess[2];
        sendExcess[0] = rank;
        sendExcess[1] = excess;
        MPI_Send(sendExcess, 2, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
        MPI_Send(slaveFlow,nodes, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);

	*/
    }
    //Master Process
    else{
        //read the input graph and generate capacity and flow matrix
        /*
        ifstream f;
        f.open("graph.txt");
        if(!f){
            cout<<"Errror opening file";
            return 0;
        }

        f >> nodes >> edges;
         */


	//Sample graph
    capacity[0][1] = 2;
    capacity[0][2] = 9;
    capacity[1][2] = 1;
    capacity[1][3] = 0;
    capacity[1][4] = 0;
    capacity[2][4] = 7;
    capacity[3][5] = 7;
    capacity[4][5] = 4;

        //flow = (double **) calloc(NODES, sizeof(double*));
        //Initialize excess map


	//printf("Capacity:\n");
    	//printMatrix(capacity);

    //printf("Max Flow:\n%d\n", pushRelabel(capacity, flow, 0, 5));

        excess[source] = INFINITE;

        //Preflow


        for(int i=0;i<NODES;i++) {		// preflow operation , (height, execess, flow)

            for(int j=0;j<NODES;j++) {

                if(i==source){
                    flow[i][j] = capacity[i][j];
                }
       	    }
        }

	//printf("Flows:\n");
   	printMatrix(flow);

        //Broadcast capacity and preflow matrix to all slaves


        MPI_Bcast(&nodes,1, MPI_INT, rank, MPI_COMM_WORLD);

        //MPI_Bcast(&capacitystring[0], (int)capacitystring.size(), MPI_CHAR, rank, MPI_COMM_WORLD);
        MPI_Bcast(&(capacity[0][0]),nodes*nodes, MPI_CHAR, rank, MPI_COMM_WORLD);
        MPI_Bcast(&(flow[0][0]), nodes*nodes, MPI_CHAR, rank, MPI_COMM_WORLD);

	cout << "Master: Broadcasted Capacity and Flow matrices \n"; // << flow[0][0];

	for(int i=1;i<nodes;i++){
		
		if(flow[0][i]!=0){		
			double t13 = flow[0][i];
			//cout << t13 << "\t";			
			MPI_Send(&t13, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
		}
	}

	cout << "Master: Sent initial flow\n";

        int got = 0;
        while(got<nodes){

            double receivedExcess[2];
            MPI_Recv(receivedExcess, 2, MPI_DOUBLE, MPI_ANY_SOURCE, 4 ,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            excess[(int)receivedExcess[0]] = receivedExcess[1];
            double receivedFlow[nodes];
            MPI_Recv(receivedFlow,nodes, MPI_DOUBLE, MPI_ANY_SOURCE, 5, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            flow[(int)receivedFlow[0]] = receivedFlow; // copy or reference copy ?
            got++;

        }

        //Discharge:
        int maxflow = 0;
        for(int i=0; i< nodes; i++)
        {
            for(int j=0; j< nodes; j++)
            {
                if(capacity[i][j]!=0)
                {
                    double delta = MIN((capacity[i][j] - flow[i][j]), excess[i]);
                    flow[i][j] += delta;
                    flow[j][i] -= delta;
                    excess[i] -= delta;
                    excess[j] += delta;
                }
            }

            if(i==sink){

                for(int k = 0;k< nodes; k++)
                {
                    if(flow[i][k]!=0)
                        maxflow += flow[i][k];
                }
            }
        }

        cout << "Maxflow is  " << maxflow;
    }





    MPI_Finalize();
    return 0;
}
