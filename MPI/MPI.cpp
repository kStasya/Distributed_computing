#include <mpi.h>
#include <cmath>


void task1() {
	MPI_Status status;
	int size, rank, recvRank;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0)
	{
		printf("Number of processes: %d.\n", size);
		printf("Hello from process %d.\n", rank);
		for (int i=1; i<size; i++)
		{
			MPI_Recv(&recvRank, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			printf("Hello from process %d.\n", recvRank);
		}
	}
	else
		MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

	MPI_Finalize();
}


void task2() {
	int rank, size, recvRank, a, b;
	double start_time, end_time;
	a = 0; b = 0;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (size != 2) {
		printf("This program requires exactly 2 MPI processes.\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	start_time = MPI_Wtime();

	if (rank==0) {
		b = 2;
		MPI_Send(&b, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
		MPI_Recv(&a, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, &status);
	}
	
	if (rank==1) {
		a = 1;
		MPI_Recv(&b, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Send(&a, 1, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
	}
	
	printf("Process %d a = %d, b = %d.\n", rank, a, b);
	end_time = MPI_Wtime();
	MPI_Finalize();
	if (rank == 0) {
		printf("Elapsed time: %fms.\n", (end_time - start_time)/10);
	}
}
	

void task3(int schema) {
	int size, rank, recvRank;
	double start_time, end_time;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// кольцо - "эстафетная палочка"
	if (schema==1) {
		if (rank == 0) {
			MPI_Send(&rank, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
			printf("Process %d sent message '%d'.\n", rank, rank);
			MPI_Recv(&recvRank, 1, MPI_INT, size-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("Process %d received message '%d'.\n", rank, recvRank);
		}
		else {
			MPI_Recv(&recvRank, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("Process %d received message '%d'.\n", rank, recvRank);
			rank = recvRank + 1;
			if (rank==size-1) {
				MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
				printf("Process %d sent message '%d'.\n", rank, rank);
			}
			else {
				MPI_Send(&rank, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
				printf("Process %d sent message '%d'.\n", rank, rank);
			}
		}
	}

	// колько - "сдвиг"
	if (schema==2) {
		MPI_Request request;

		int prev = rank - 1; 
		int next = rank + 1;
		if (rank==0) 
			prev = size - 1;
		if (rank==size-1)
			next = 0;

		MPI_Isend(&rank, 1, MPI_INT, next, 0, MPI_COMM_WORLD, &request);
		MPI_Irecv(&recvRank, 1, MPI_INT, prev, 0, MPI_COMM_WORLD, &request);
		MPI_Wait(&request, &status);
		for (int i=0; i<size; i++) {
			if (rank == i) {
				printf("Process %d sent message '%d' to %d.\n", rank, rank, next);
				printf("Process %d received message '%d' from %d.\n", rank, recvRank, prev);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

	// master-slave
	if (schema==3) {
		if (rank==0)
		{
			for (int i=1; i<size; i++)
			{
				MPI_Recv(&recvRank, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				printf("Process %d received a message '%d' from %d.\n", rank, recvRank, recvRank);
			}
		}
		else
		MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
    
	// каждый -> каждому
	if (schema==4) {
		int send[size];
		int recv[size];
		int i;
		for (i=0; i<size; i++) {
			send[i] = rank;
		}
		for (i=0; i<size; i++) {
			MPI_Send(&send[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		}
		for (i=0; i < size; i++) {
			MPI_Recv(&recv[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("Process %d received message %d from process %d.\n", rank, recv[i], i);
		}
	}
	MPI_Finalize();
}


long double maintask(int rank, int size, int N, float R) {
	long double answer, x, y, z, loc_sum, sum=0, V_sphere = 4 / 3 * M_PI * R * R;

	// calculating volume of one octant
	int n = 0;
	for(int i=rank; i<=N; i+=size){
		for(int j=0; j<=N; j++){
			for(int k=0; k<=N; k++){
				x = R * i / N;
				y = R * j / N;
				z = R * k / N;
				// checking inner points
				if (x*x+y*y+z*z<=R*R)
				{
					n++;
					sum += exp(-x*x - y*y - z*z);
				}
			}
		}   
	}
	MPI_Reduce(&sum, &loc_sum, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	answer = 8 * V_sphere * sum / n;
	return answer;
}

int main (int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	//task1();
	//task2();
	//task3(1);
	//task3(2);
	//task3(3);
	//task3(4);

	/* main task */
	float eps = 1e-5, R = 1;
	int rank, size, N = 500;
	double start_time, end_time;
	long double I, newI = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	printf('Number of processes: %d.\n', size);
	start_time = MPI_Wtime();
	do {
		I = maintask(rank, size, N, R);
		N *= 2;
		newI = maintask(rank, size, N, R);
		if (rank == 0) {
			printf("Estimated error: %f.\n", fabs(I - newI));
			if (fabs(I - newI)>=eps) {
				printf("Integral: %Lf, N: %d.\n", newI, N);
			}
			else{
				end_time = MPI_Wtime();
				printf("Integral: %Lf, N: %d, elapsed time: %fs\n", newI, N, end_time - start_time);
				MPI_Finalize();
				return 0;
			}
		}
	} while (true);
}
