#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N 1<<3

void allocate(
    double **vector_pointer,
    int size);

void fill(
    double vector_pointer[],
    size_t size);

void print(
    double *vector_pointer,
    size_t size);

void sum(
    double *vector_pointer1,
    double *vector_pointer2,
    double *vector_pointer3,
    size_t size);

int main(int argc, char *argv[]){        

    double *x = NULL, *y = NULL;
    double *local_x = NULL, *local_y = NULL;

    int comm_sz;
    int my_rank;
    int local_n;

    int n = N;
    int f = 0;

    if (argc > 1){
        n = 1 << atoi(argv[1]);
        if (argc > 2) {
            f = atoi(argv[2]);
        }
    }

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    local_n = n / comm_sz;

    allocate(&local_x, local_n);
    allocate(&local_y, local_n);

    double starttime, endtime;
       
    if (my_rank == 0){
        
        allocate(&x, n);
        allocate(&y, n);
        fill(x, n);
        fill(y, n);

        starttime = MPI_Wtime();

    }
      
    double scatter_time = MPI_Wtime();
    MPI_Scatter(x, local_n, MPI_DOUBLE, local_x, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(y, local_n, MPI_DOUBLE, local_y, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double scatter_end_time = MPI_Wtime();
    sum(local_x, local_y, local_y, local_n);
    double gather_time = MPI_Wtime();
    MPI_Gather(
        local_y, 
        local_n, 
        MPI_DOUBLE, 
        y, 
        local_n, 
        MPI_DOUBLE, 
        0, 
        MPI_COMM_WORLD
    );
    double gather_end_time = MPI_Wtime();
    
    if (my_rank == 0){

        if(f == 0) print(y, n);
        
        free(x);
        free(y);

        endtime   = MPI_Wtime();

        printf("That took %f seconds\n",endtime-starttime);
        printf("Tempo gasto para o scatter: %lf\n", scatter_end_time - scatter_time);
        printf("Tempo gasto para o gather: %lf\n", gather_end_time - gather_time);
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
    size_t size){
    for (size_t i = 0; i < size; i++)
        vector_pointer[i] = (1 + i)/2;
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
