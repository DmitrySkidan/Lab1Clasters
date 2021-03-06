#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
const int MY_TAG = 42;
int main(int argc, char* argv[]){
    int rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0) {
        sleep(5);
        int data[] = { 1, 2, 3, 4, 5};
        MPI_Request send_req;
        MPI_Isend(&data, 5, MPI_INT, 1, MY_TAG, MPI_COMM_WORLD, &send_req);
        sleep(5);
        MPI_Wait(&send_req, MPI_STATUS_IGNORE);
    }
    if(rank == 1){/*Процес 1 починає отримувати дані раніше ніж процес 0 почав їх передавати*/
        int data[5];
        MPI_Request recv_req;
        MPI_Irecv(&data, 5, MPI_INT, 0, MY_TAG, MPI_COMM_WORLD, &recv_req); /* буфер data недоступний… */
        sleep(2);
        MPI_Wait(&recv_req, MPI_STATUS_IGNORE);
        for(int i = 0; i < 5; i++) {
            printf(" %d", data[i]);
        }
    }
    MPI_Finalize();
    return 0;
}
