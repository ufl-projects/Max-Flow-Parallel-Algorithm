#include <fstream>

#include <stdlib.h>

#include <stdio.h>

#include <sstream>

#include <sys/time.h>

#include "mpi.h"



#define NODES 8

#define SINK 7

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

#define INFINITE 10000000

#define PARALLEL


#define DEBUG


using namespace std;





#ifdef PARALLEL

void push(const double * const * C, double ** F, double *excess, int u, int v) {

    int send = MIN(excess[u], C[u][v] - F[u][v]);

    F[u][v] += send;

    F[v][u] -= send;

    excess[u] -= send;

    excess[v] += send;

}


void relabel(const double * const * C, const double * const * F, int *height, int u) {

    int v;

    int min_height = INFINITE;

    for (v = 0; v < NODES; v++) {

        if (C[u][v] - F[u][v] > 0) {

            min_height = MIN(min_height, height[v]);

            height[u] = min_height + 1;

        }

    }

};


void discharge(const double * const * C, double ** F, double *excess, int *height,double *seen, int u) {

    while (excess[u] > 0) {

        if (seen[u] < NODES) {

            int v = seen[u];

            if ((C[u][v] - F[u][v] > 0)&& (height[u] > height[v])){

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


void printMatrix(const double * const * M) {

    int i,j;

    for (i = 0; i < NODES; i++) {

        for (j = 0; j < NODES; j++)

            printf("%f\t",M[i][j]);

        printf("\n");

    }

}



double **alloc_2d_int(int rows, int cols) {

    double *data = (double *)calloc(rows*cols, sizeof(double));

    double **array= (double **)calloc(rows,sizeof(double*));

    for (int i=0; i<rows; i++)

        array[i] = &(data[cols*i]);

    

    return array;

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

    

    char processor_name[MPI_MAX_PROCESSOR_NAME];

    int name_len;

    MPI_Get_processor_name(processor_name, &name_len);

    

    string capacitystring;

    string flowstring;

    

    

    int i,nodes,edges;

    int source = 0;

    int sink=SINK;

    nodes=NODES;

    

#ifdef DEBUG

    cout<<"Initialized MPI "<<endl;

#endif

    

    

    

    

#ifdef DEBUG

    cout<<"Master/Slave operation begins "<<endl;

#endif

    

    

    //Slave process

    if(rank!=0){

        

       

        

        

#ifdef DEBUG

        cout<<"Slave "<<rank<<" : is receiving the broadcast info"<<endl;

#endif

        

        

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Bcast(&nodes,1, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        //MPI_Recv(&node,1,MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        double **flow, **capacity;

        capacity = alloc_2d_int(nodes, nodes);

        flow = alloc_2d_int(nodes, nodes);

        double excess=0,push_flow=0;

#ifdef DEBUG

        cout<<"Slave "<<rank<<" : No of nodes received -  "<< nodes<<endl;

#endif

        //cout <<  "My nodes " << nodes << "\n";

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Bcast(&(capacity[0][0]),nodes*nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Bcast(&(flow[0][0]), nodes*nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

	
        
        double* slaveFlow = new double[nodes+1]();

        

        

        //Receive excess from adjacent nodes

#ifdef DEBUG

        cout<<"Slave "<<rank<<" : Receive excess from adjacent nodes"<<endl;

#endif

        for(int i = 0; i<rank; i++)

        {

            if( capacity[i][rank] != 0 ) {

                double t;

                MPI_Recv(&t,1,MPI_DOUBLE, i, 2, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

                slaveFlow[i] = -t;

                excess += t;

#ifdef DEBUG

                cout<<"Slave "<<rank<<" : EXCESS  value: "<<excess<<endl;

#endif

                

            }

        }

        

#ifdef DEBUG

        cout<<"Slave "<<rank<<" : Calculating excess and pushing forward"<<endl;

#endif

        for(int i = rank+1; i<=sink;i++)

        {

            if ( capacity[rank][i] != 0 && excess != 0 )

            {

                push_flow = MIN( capacity[rank][i] , excess );

                slaveFlow[i] = push_flow;

                excess -= push_flow;

                

                MPI_Send(&push_flow, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);

            }

            

            if(excess==0) MPI_Send(&excess, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);

        }

        

        

#ifdef DEBUG

        cout<<"Slave "<<rank<<" : Calculating excess and pushing forward"<<endl;

#endif
	
        double sendExcess[2];

        sendExcess[0] = rank;

        sendExcess[1] = excess;

        MPI_Send(&sendExcess[0], 2, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);

        MPI_Send(&slaveFlow[0],nodes, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);

#ifdef DEBUG

        cout<<"Slave "<<rank<<" : **Slave flow matrix***"<<endl;

#endif
	
        cout<<"Slave "<<rank<<" : **Slave flow matrix***"<<endl;

        //Print slaveFlow

        

         for(int j = 0; j < NODES; j++) {

         	printf("\n");

         	printf("%f ", slaveFlow[j]);

         }

    }

    //Master Process

    else{

        

        struct timeval start, end;

        gettimeofday(&start, NULL);

        

        // benchmark code

        

        //read the input graph and generate capacity and flow matrix     

        

        //flow = (double **) calloc(NODES, sizeof(double*));

        

        //Initialize excess map

        

        

#ifdef DEBUG

        cout<<"Master initiates broadcasting"<<endl;

#endif

        double **flow, **capacity,*excess,*seen;

        excess = (double*) calloc(NODES, sizeof(double));

        seen = (double*) calloc(NODES, sizeof(double));

        

        capacity = alloc_2d_int(nodes, nodes);

        flow = alloc_2d_int(nodes, nodes);

        

        

        

        

        int *height = (int *) calloc(NODES, sizeof(int));

        height[source] = NODES;

        int *list = (int *) calloc((NODES-2), sizeof(int));

        

        

        double capacityValue;

        int n1,n2;

        ifstream f;	//("t1small.txt");	//"facebook_combined.txt");

        f.open("small-graph.txt");

        if(!f){

            cout<<"Errror opening file";

            return 0;

        }

	cout << "File opened \n";

        f >> nodes;

	cout << "Nodes " << nodes;

        f >> edges;

        f >> n1 >> n2 >> capacityValue;

        while(!f.eof())

        {

            capacity[n1][n2]=capacityValue;

            f >> n1 >> n2 >> capacityValue;

        }

	cout << "File read \n";

        f.close();


        for (int i = 0, p = 0; i < NODES; i++){

            if((i != source) && (i != sink)) {

                list[p] = i;

                p++;

            }

        }

        //excess[source] = INFINITE;

        

#ifdef DEBUG

        cout<<"Master:Initializing capacity and flow matrix"<<endl;

#endif

        

        

        //Sample graph

        /*

         capacity[0][1] = 1;

         capacity[0][2] = 2;

         capacity[1][2] = 1;

         capacity[1][3] = 0;

         capacity[1][4] = 0;

         capacity[2][4] = 7;

         capacity[3][5] = 7;

         capacity[4][5] = 4;

         

        capacity[0][1] = 3;

        capacity[0][3] = 3;

        capacity[1][2] = 4;

        capacity[2][0] = 3;

        capacity[2][3] = 1;

        capacity[2][4] = 2;

        capacity[3][4] = 2;

        capacity[3][5] = 6;

        capacity[4][1] = 1;

        capacity[4][6] = 1;

        capacity[5][6] = 9;

        

         capacity[1][2] = 3;

         capacity[1][3] = 3;

         capacity[2][4] = 1;

         capacity[2][5] = 1;

         capacity[3][4] = 2;

         capacity[4][5] = 2;

         

         

         capacity[0][1] = 1;

         capacity[0][2] = 2;

         capacity[1][2] = 1;

         capacity[1][3] = 0;

         capacity[1][4] = 0;

         capacity[2][4] = 7;

         capacity[3][5] = 7;

         capacity[4][5] = 4;

         */

        printf("Capacity:\n");

        printMatrix(capacity);

        //printMatrix(flow);

        //printf("Max Flow:\n%d\n", pushRelabel(capacity, flow, 0, 5));

        

        

        

        

        

#ifdef DEBUG

        cout<<"Master : "<<rank<<" Initiated Preflow operation on Master"<<endl;

#endif

        

        

        //Preflow

        

        

        for(int i=0;i<NODES;i++){

            //push(capacity, flow, excess, source, i);

            for(int j=0;j<NODES;j++) {

                

                if(i==source){

                    flow[i][j] = capacity[i][j];

                }

            }

        }

        

        

        //Broadcast capacity and preflow matrix to all slaves

        

        printf("Flows:\n");

        printMatrix(flow);

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Bcast(&nodes, 1, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Bcast(&(capacity[0][0]),nodes*nodes, MPI_DOUBLE, rank, MPI_COMM_WORLD);

        MPI_Bcast(&(flow[0][0]), nodes*nodes, MPI_CHAR, rank, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        

        

#ifdef DEBUG

        cout<<"Master : "<<rank<<" Sending flow info for source nodes"<<endl;

#endif

        

        

        //Send flow info

        for(int i=1;i<nodes;i++){

            

            if(flow[0][i]!=0){

                double t13 = flow[0][i];

                //cout << t13 << "\t";

                MPI_Send(&t13, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);

            }

        }

        

        

#ifdef DEBUG

        cout<<"Master : "<<rank<<" Sent flow info for source nodes"<<endl;

#endif

        

        

        int got = 1;

        while(got<nodes){

		cout << "start \n";
	
            double receivedExcess[2];

            

            MPI_Recv(receivedExcess, 2, MPI_DOUBLE, got, 4 ,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

            int ndex = (int)receivedExcess[0];

            double excess12 = receivedExcess[1];

            excess[ndex] = excess12;

            

            

#ifdef DEBUG

            cout<<"Master : "<<rank<<" Excess of "<<ndex<<" is  "<<excess12<<endl;

#endif
	   double receivedFlow[nodes];

            MPI_Recv(receivedFlow,nodes, MPI_DOUBLE, got, 5, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

            //flow[ndex] = receivedFlow;

            for(int i=0; i< nodes; i++)

            { flow[ndex][i] = receivedFlow[i]; }

#ifdef DEBUG

            cout<<"Master : "<<rank<<" Flow of "<<got  <<endl;

            //printMatrix(flow);

#endif

	
            got++;

           cout << "end \n";

        }

        
        

         int p = 0;

         while (p < NODES - 2) {

		cout << "DOne";

		 int u = list[p];

		 int old_height = height[u];

		 discharge(capacity, flow, excess,height, seen, u);

		 if (height[u] > old_height) {

			 moveToFront(p,list);

			 p=0;

		 }

		 else

		 	p += 1;

         }

        
        int maxflow = 0;

        for (i = 0; i < NODES; i++)

            maxflow += flow[source][i];

        cout << "Maxflow is  " << maxflow<<endl;

        

        for(int j = 0; j < NODES; j++) {

            printf("%f ", excess[j]);

        }

        printf("\n");

        /*

         int maxflow = 0;

         for (i = 0; i < NODES; i++)

         maxflow += flow[source][i];

         cout << "Maxflow is  " << maxflow<<endl;

         for(int j = 0; j < NODES; j++) {

         printf("%f ", excess[j]);

         }*/

        printf("\n");

        //cout<<"Flow";

        printMatrix(flow);

        

        

        gettimeofday(&end, NULL);

        

        double executionTime = ((end.tv_sec  - start.tv_sec) * 1000000u +

                                end.tv_usec - start.tv_usec) / 1.e6;

        

        cout<<"EXECUTION TIME : "<<executionTime<<" secs"<<endl;


        

    }

    

    

    

    MPI_Finalize();

    return 0;

}

#endif
