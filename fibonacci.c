//
// Created by Zhe Wang on 11/22/2024.
//

#include <stdio.h>

#include "mpi.h"

void fibonacci(int my_id, int size) {
    int num1 = my_id;
    int num2 = 0;
    for (int i = my_id - 1; i >= 0 && i >= my_id - 2; i--) {
        MPI_Recv(i == my_id - 2 ? &num1 : &num2, 1, MPI_INT, i, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    const int sum = num1 + num2;
    printf("id %d: %d\n", my_id, sum);
    for (int i = my_id + 1; i < size && i <= my_id + 2; i++) {
        MPI_Send(&sum, 1, MPI_INT, i, 99, MPI_COMM_WORLD);
    }
}

int main(int argc, char *argv[]) {
    int my_id = 0;
    int size = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    fibonacci(my_id, size);

    if (my_id == size - 1) {
        putchar('\n');
    }
    MPI_Finalize();
    return 0;
}