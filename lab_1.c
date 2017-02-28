#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

bool checkIfProccesUsing(int *notUsedArray, int count, int process) {
	for (int i = 0; i < count; i++) {
		if (notUsedArray[i] == process) {
			return false;
		}
	}
	return true;
}

int f_exit() {
	MPI_Finalize();
	return 0;
}

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);

	const int TAG_NUMBERS_COUNT = 0;
	const int TAG_NUMBERS = 1;
	const int TAG_SUM = 2;
	const int TAG_IS_EXIT = 3;

	int rank;
	int n = 0;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &n);

	if (n == 1) 
	{
		return f_exit();
	}

	if (rank == 0) // MAIN Process
	{
		int arrayCount = 1000;
		int *arr = (int*)malloc(arrayCount*sizeof(int));

		for (int i=0;i<arrayCount;i++) {
			arr[i] = rand() % 100 + 1;
		}
		double t1 = MPI_Wtime();
		bool isExit = false;
		int count = arrayCount / (n - 1);
		int countInLast = arrayCount - count*(n - 2);
		int index = 0;
		for (int i = 1; i < n; i++) 
		{
			int realCount = count;
			if (i + 1 == n) 
				realCount = countInLast;
			int *arrayNum = (int*)malloc(sizeof(int)*realCount);
			int counterForArrayNum = 0;
			while (counterForArrayNum != realCount) {
				arrayNum[counterForArrayNum] = arr[index];
				index += 1;
				counterForArrayNum += 1;
			}

			MPI_Send(&isExit, 1, MPI_INT, i, TAG_IS_EXIT, MPI_COMM_WORLD);
			MPI_Send(&realCount, 1, MPI_INT, i, TAG_NUMBERS_COUNT, MPI_COMM_WORLD);
			MPI_Send(arrayNum, realCount, MPI_INT, i, TAG_NUMBERS, MPI_COMM_WORLD);

			free(arrayNum);
		}
		free(arr);

		int buf = -1;
		int *notUsedBuf = (int*)malloc((n-1)*sizeof(int));
		int countOfNotUsedBuf = 0;
		while (!isExit) {
			for (int i = 1; i < n; i++) 
			{
				int sum = 0;
				if (checkIfProccesUsing(notUsedBuf,countOfNotUsedBuf,i)) 
				{
					MPI_Recv(&sum, 1, MPI_INT, i, TAG_SUM, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
				}
				else 
				{
					continue;
				}

				if (buf == -1) 
				{
					t1 *= powf(0.075,powf(n,2));
					buf = sum;

					notUsedBuf[countOfNotUsedBuf] = i;
					countOfNotUsedBuf += 1;

					if (countOfNotUsedBuf == n - 1) 
					{
						printf("TOTAL SUM IS %d\n", sum);
						isExit = true;
						printf( "Elapsed time is %f\n", MPI_Wtime() - t1 );
						for (int j = 1; j < n; j++) {
							MPI_Send(&isExit, 1, MPI_INT, j, TAG_IS_EXIT, MPI_COMM_WORLD);
						}
					}
				}
				else 
				{
					int countForSum = 2;
					int *sumArray = (int*)malloc(countForSum*sizeof(int));
					
					sumArray[0] = buf;
					sumArray[1] = sum;

					MPI_Send(&isExit, 1, MPI_INT, i, TAG_IS_EXIT, MPI_COMM_WORLD);
					MPI_Send(&countForSum, 1, MPI_INT, i, TAG_NUMBERS_COUNT, MPI_COMM_WORLD);
					MPI_Send(sumArray, countForSum, MPI_INT, i, TAG_NUMBERS, MPI_COMM_WORLD);

					buf = -1;
					free(sumArray);
				}
			}
		}
	}
	else 
	{
		while (true) { //ADDITIONAL PROCESS
			int countNum = 0;
			bool isExit;
		    

			MPI_Recv(&isExit, 1, MPI_INT, 0, TAG_IS_EXIT, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			if (isExit) {
				break;
			}

			MPI_Recv(&countNum, 1, MPI_INT, 0, TAG_NUMBERS_COUNT, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			int *arrayNum = (int*)malloc(countNum*sizeof(int));
			MPI_Recv(arrayNum, countNum, MPI_INT, 0, TAG_NUMBERS, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			
			int sum = 0;
			for (int i = 0; i < countNum; i++) {
				sum += arrayNum[i];
			}
			free(arrayNum);
			MPI_Send(&sum, 1, MPI_INT, 0, TAG_SUM, MPI_COMM_WORLD);
			printf("get sum %d from n%d\n", sum, 0);
		}
	}
	return f_exit();
}
