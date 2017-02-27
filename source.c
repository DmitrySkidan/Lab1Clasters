#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

void addElementToArray(int **arrayOf, int n, int element) {
	int *copyarray = (int*)malloc(sizeof(int)*(n));
	for (int i = 0; i < n - 1; i++) {
		copyarray[i] = (*arrayOf)[i];
	}
	copyarray[n - 1] = element;
	free(*arrayOf);
	*arrayOf = copyarray;
}

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
	const int TAG_IS_EXIT_DONE = 4;

	int rank;
	int n = 0;

	int countOfUsage = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &n); 

	// добавить проверку, добавить генерацию, получение аргумента со строки

	if (n == 1) 
	{
		return f_exit();
	}

	if (rank == 0) // MAIN Process
	{
		int arrayCount = 100;
		int *arr = (int*)malloc(sizeof(int)*arrayCount);

		for (int i=0;i<arrayCount;i++) {
			arr[i] = rand() % 100 + 1;
		}
		clock_t start = clock(), diff;
		bool isExit = false;
		int count = arrayCount / (n - 1);
		int countInLast = arrayCount - count*(n - 2);
		int index = 0;
		for (int i = 1; i < n; i++) 
		{
			int realCount = count;
			if (i + 1 == n) 
				realCount = countInLast;
			int *arrayNum = (int*)malloc(sizeof(int));
			int counterForArrayNum = 1;
			while (counterForArrayNum != realCount + 1) {
				addElementToArray(&arrayNum, counterForArrayNum, arr[index]);
				index += 1;
				counterForArrayNum += 1;
			}

			MPI_Send(&isExit, 1, MPI_INT, i, TAG_IS_EXIT, MPI_COMM_WORLD);
			MPI_Send(&realCount, 1, MPI_INT, i, TAG_NUMBERS_COUNT, MPI_COMM_WORLD);
			MPI_Send(arrayNum, realCount, MPI_INT, i, TAG_NUMBERS, MPI_COMM_WORLD);
		}
		
		int buf = -1;
		int *notUsedBuf = (int*)malloc(sizeof(int));
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
					buf = sum;

					countOfNotUsedBuf += 1;
					addElementToArray(&notUsedBuf,countOfNotUsedBuf,i);

					if (countOfNotUsedBuf == n - 1) 
					{
						printf("TOTAL SUM IS %d\n", sum);
						isExit = true;
						int msec = diff * 1000 / CLOCKS_PER_SEC;
						printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
						for (int j = 1; j < n; j++) {
							MPI_Send(&isExit, 1, MPI_INT, j, TAG_IS_EXIT, MPI_COMM_WORLD);
						}
					}
				}
				else 
				{
					int *sumArray = (int*)malloc(sizeof(int));
					int countForSum = 2;
					addElementToArray(&sumArray,1,buf);
					addElementToArray(&sumArray,2,sum);

					MPI_Send(&isExit, 1, MPI_INT, i, TAG_IS_EXIT, MPI_COMM_WORLD);
					MPI_Send(&countForSum, 1, MPI_INT, i, TAG_NUMBERS_COUNT, MPI_COMM_WORLD);
					MPI_Send(sumArray, countForSum, MPI_INT, i, TAG_NUMBERS, MPI_COMM_WORLD);

					buf = -1;
				}
			}
		}
	}
	else 
	{
		while (true) { //ADDITIONAL PROCESS
			int countNum = 0;
			bool isExit;
		    int *arrayNum = (int*)malloc(sizeof(int));

			MPI_Recv(&isExit, 1, MPI_INT, 0, TAG_IS_EXIT, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			if (isExit) {
				break;
			}

			MPI_Recv(&countNum, 1, MPI_INT, 0, TAG_NUMBERS_COUNT, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(arrayNum, countNum, MPI_INT, 0, TAG_NUMBERS, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			
			int sum = 0;
			for (int i = 0; i < countNum; i++) {
				sum += arrayNum[i];
			}
			
			MPI_Send(&sum, 1, MPI_INT, 0, TAG_SUM, MPI_COMM_WORLD);
			printf("get sum %d from n%d\n", sum, 0);
		}
	}
	return f_exit();
}