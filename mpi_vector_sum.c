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

int main(int argc, char** argv){        

    double *x = NULL, *y = NULL;
    double *local_x = NULL, *local_y = NULL;
    long vec_size = strtol(argv[1], NULL, 10);
    int gCount = 0;

    int comm_sz;
    int my_rank;
    int local_n;
    int cond = 0;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    /*if(vec_size % comm_sz == 0) {
        cond = 1;
        local_n = vec_size / comm_sz;
    }
    else {
        local_n = (int)vec_size / comm_sz;
        if(my_rank == 0) local_n = vec_size % comm_sz;
    }*/
    local_n = vec_size / comm_sz;
    allocate(&local_x, local_n);
    allocate(&local_y, local_n);

    double starttime, endtime;
    if (my_rank == 0){   
        allocate(&x, vec_size);
        allocate(&y, vec_size);
        fill(x, 0.1, vec_size);
        fill(y, 0.1, vec_size);
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
