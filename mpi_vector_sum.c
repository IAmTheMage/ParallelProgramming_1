/* File:     vector_sum.c
 *
 * Purpose:  Implement vector addition
 *
 * Compile:  gcc -g -Wall -o vector_add vector_sum.c
 * Run:      ./vector_sum
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N 1<<3

void allocate(
    double **vector_pointer,
    int size);

void fill(
    double vector_pointer[],
    double ratio,
    size_t size);

void print(
    double *vector_pointer,
    size_t size);

void sum(
    double *vector_pointer1,
    double *vector_pointer2,
    double *vector_pointer3,
    size_t size);

void allocate_disp(int** disp, size_t size);
void allocate_send_counts(int** sendCounts, size_t size);

int main(int argc, char** argv){        

    double *x = NULL, *y = NULL;
    double *local_x = NULL, *local_y = NULL;
    double* scatterV_Test;
    int local_received = 0;
    int mod = 0;
    long vec_size = strtol(argv[1], NULL, 10);
    int* disp;
    int* sendCounts;

    int comm_sz;
    int my_rank;
    int local_n;
    int cond = 0;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    local_n = (int)vec_size / comm_sz;
    local_received = local_n;
    mod = vec_size % comm_sz;
    if(my_rank == 0) {
        local_received += mod;
    }
    allocate(&local_x, local_received);
    allocate(&local_y, local_received);
    double starttime, endtime;
    if (my_rank == 0){   
        allocate(&x, vec_size);
        allocate(&y, vec_size);
        fill(x, 1.0, vec_size);
        fill(y, 1.0, vec_size);
        starttime = MPI_Wtime();

    }
      
    allocate_disp(&disp, comm_sz);
    allocate_send_counts(&sendCounts, comm_sz);
    for(int i = 0; i < comm_sz; i++) {
        if(i == 0) {
            sendCounts[i] = local_n + mod;
            disp[i] = 0;
        }
        else {
            sendCounts[i] = local_n;
            disp[i] = i * local_n + mod;
        }
    }
    double scatterv_time = MPI_Wtime();
    MPI_Scatterv( 
        x, 
        sendCounts, 
        disp,
        MPI_DOUBLE, 
        local_x, 
        local_received, 
        MPI_DOUBLE,
        0, 
        MPI_COMM_WORLD
    );
    MPI_Scatterv( 
        y, 
        sendCounts, 
        disp,
        MPI_DOUBLE, 
        local_y, 
        local_received, 
        MPI_DOUBLE,
        0, 
        MPI_COMM_WORLD
    );
    double scatterv_end_time = MPI_Wtime();
    sum(local_x, local_y, local_y, local_received);
    double gatherv_time = MPI_Wtime();
    MPI_Gatherv( 
        local_y, 
        local_received, 
        MPI_DOUBLE, 
        y, 
        sendCounts, 
        disp, 
        MPI_DOUBLE, 
        0, MPI_COMM_WORLD
    );
    double gatherv_end_time = MPI_Wtime();
    if (my_rank == 0){
        // print(y, N);
        free(x);
        free(y);

        endtime   = MPI_Wtime();
        printf("That took %f seconds\n",endtime-starttime);
        printf("Tempo gasto para o scatterv: %lf\n", scatterv_end_time - scatterv_time);
        printf("Tempo gasto para o gatherv: %lf\n", gatherv_end_time - gatherv_time);
    }

    free(local_x);
    free(local_y);
    
    MPI_Finalize();

    return 0;
}

void allocate(
    double **vector_pointer,
    int size){
    *vector_pointer = (double *)malloc(size * sizeof(double));
    if (vector_pointer == NULL)
    {
        fprintf(stderr, "Can't allocate vector\n");
        exit(-1);
    }
}

void fill(
    double vector_pointer[],
    double ratio,
    size_t size){
    for (size_t i = 0; i < size; i++)
        vector_pointer[i] = ratio * i;
}


void allocate_send_counts(
    int** sendCounts, size_t size
)
{
    *sendCounts = (int*)malloc(size * sizeof(int));
}

void allocate_disp(int** disp, size_t size) {
    *disp = (int*)malloc(size * sizeof(int));
}

void print(
    double *vector_pointer,
    size_t size){
    for (size_t i = 0; i < size; i++)
        printf("%f ", vector_pointer[i]);
}

void sum(
    double *vector_pointer1,
    double *vector_pointer2,
    double *vector_pointer3,
    size_t size){
    for (size_t i = 0; i < size; i++)
        vector_pointer3[i] = vector_pointer1[i] + vector_pointer2[i];
}
