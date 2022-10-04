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
#include </home/ronaldo/opt/openmpi/include/mpi.h>

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

int main(void){        

    double *x = NULL, *y = NULL;
    double *local_x = NULL, *local_y = NULL;

    int comm_sz;
    int my_rank;
    int local_n;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    local_n = N / comm_sz;

    allocate(&local_x, local_n);
    allocate(&local_y, local_n);

    double starttime, endtime;
       
    if (my_rank == 0){
        
        allocate(&x, N);
        allocate(&y, N);
        fill(x, 0.1, N);
        fill(y, 0.1, N);

        starttime = MPI_Wtime();

    }
      

    MPI_Scatter(x, local_n, MPI_DOUBLE, local_x, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(y, local_n, MPI_DOUBLE, local_y, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    sum(local_x, local_y, local_y, local_n);

    MPI_Gather(local_y, local_n, MPI_DOUBLE, y, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (my_rank == 0){

        // print(y, N);
        
        free(x);
        free(y);

        endtime   = MPI_Wtime();
        printf("That took %f seconds\n",endtime-starttime);

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
