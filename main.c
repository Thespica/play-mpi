/*
 * Copyright 2024 Thespica
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#include "mpi.h"

int sum_to_num(int sum_to) {
    int my_id = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int need_sum_to = 100;
    int sum = 0;

    const int data_per_process = (need_sum_to + size - 1) / size;
    const int begin = my_id * data_per_process + 1;
    int end = begin + data_per_process;
    for (int i = begin; i < end && i <= need_sum_to; i++) {
        sum += i;
    }
    printf("pid: %d: %d\n", my_id, sum);
    if (my_id == 0) {
        int receive = 0;
        for (int recv_pid = 1; recv_pid < size; recv_pid++) {
            MPI_Recv(&receive, 1, MPI_INT, recv_pid, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum += receive;
        }
        printf("all sum: %d\n", sum);
        return sum;
    }
    if (my_id == 1) {
        MPI_Send(&sum, 1, MPI_INT, 0, 99, MPI_COMM_WORLD);
    }
    return 0;
}

void block(void) {
    int my_id = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    srand(time(NULL) + my_id);
    int rand_value = rand();
    if (my_id == 0) {
        int target_pid = 1;
        int other_random = 0;
        MPI_Send(&rand_value, 1, MPI_INT, target_pid, 99, MPI_COMM_WORLD);
        MPI_Recv(&other_random, 1, MPI_INT, target_pid, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int min_value = rand_value < other_random ? rand_value : other_random;
        printf("min random value from id %d: %d\n", my_id, min_value);
    }
    if (my_id == 1) {
        int target_pid = 0;
        int other_random = 0;
        MPI_Send(&rand_value, 1, MPI_INT, target_pid, 99, MPI_COMM_WORLD);
        MPI_Recv(&other_random, 1, MPI_INT, target_pid, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int max_value = rand_value > other_random ? rand_value : other_random;
        printf("max random value from id %d: %d\n", my_id, max_value);
    }
}

void hava_dependency() {
    int my_id = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int sum = my_id == 0 ? 10 : 1;
    if (my_id > 0) {
        int receive = 0;
        MPI_Recv(&receive, 1, MPI_INT, my_id - 1, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        sum += receive;
    }
    if (my_id < size - 1) {
        MPI_Send(&sum, 1, MPI_INT, my_id + 1, 99, MPI_COMM_WORLD);
    }
    printf("id %d sum: %d\n", my_id, sum);
}

void my_broadcast(double *buff, int count, int root) {
//    srand(time(NULL));
//    double rds[5];
//    if (my_id == 0) {
//        for (int i = 0; i < 5; i++) {
//            rds[i] = rand() / (RAND_MAX + 1.0);
//        }
//    }
//    my_broadcast(rds, 5, 0);

    int my_id = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (my_id == 0) {
        for (int i = 1; i < size; i++) {
            MPI_Send(buff, count, MPI_DOUBLE, i, 99, MPI_COMM_WORLD);
        }
    }
    if (my_id != 0) {
        MPI_Recv(buff, count, MPI_DOUBLE, root, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    printf("my_id=%d, rd=:\n", my_id);
    for (int i = 0; i < count; i++) {
        printf("%f ", buff[i]);
    }
    putchar('\n');
}

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

void calc_pi(int my_id, int size) {
    double sum = 0;
#define N 10000
    for (int i = my_id; i < N; i += size) {
        sum += (1 / N) * (4 / (1 + pow(i / N, 2)));
    }
    double global_sum = 0;
    printf("%d proc: %f\n", my_id, sum);
    MPI_Reduce(&sum,       // 输入缓冲区（本地值）
               &global_sum,        // 输出缓冲区（全局总和，只有 root 进程有值）
               1,                  // 输入缓冲区中元素个数
               MPI_DOUBLE,            // 数据类型
               MPI_SUM,            // 操作类型（求和）
               0,                  // root 进程的 rank
               MPI_COMM_WORLD);    // 通信域
    if (my_id == 0) {
        double pi = 4 * (1 - sum);
        printf("sum: %f\n", pi);
    }
}

void monte_carlo_pi(int my_id, int size) {
    srand(time(NULL) + my_id);
    int m = 100000;
    int n = 0;
    for (int i = 0; i < m; i++) {
        double x = (double)rand() / (RAND_MAX + 1.0);
        double y = (double)rand() / (RAND_MAX + 1.0);
        n += (pow(x - 0.5, 2) + pow(y - 0.5, 2) < 0.25);
    }
    printf("%d proc: %f\n", my_id, 4.0F * n / m);
    int global_m;
    int global_n;
    MPI_Reduce(&m, &global_m, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&n, &global_n, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (my_id == 0) {
        printf("all sum: %f\n", 4.0F * global_n / global_m);
//        printf("all m: %d", global_m);
//        printf("all n: %d", global_n);
    }
}

void calc_pi_by_expression(int my_id, int size, int n) {
    double sum_of_single = 0;
    for (int i = my_id; i <= n; i += size) {
        sum_of_single += (1.0 / pow(16, i)) * ((4.0 / (8 * i + 1)) -
                                               (2.0 / (8 * i + 4)) -
                                               (1.0 / (8 * i + 5)) -
                                               (1.0 / (8 * i + 6)));
    }
    double receive_sum = 0;
    MPI_Reduce(&sum_of_single, &receive_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (my_id == 0) {
        printf("all sum: %f\n", receive_sum);
    }
}

void calc_e(int my_id, int size, int n) {
    double single_sum = 0;
    double reciprocal_of_factorial = 1;
    for (int i = 2; i <= my_id; i++) {
        reciprocal_of_factorial *= 1.0 / i;
    }
    for (int i = my_id; i <= n; i += size) {
        for (int j = i - size + 1; j <= i && j >= 1; j++) {
            reciprocal_of_factorial *= 1.0 / j;
        }
        single_sum += reciprocal_of_factorial;
    }
    double all_sum = 0;
    MPI_Reduce(&single_sum, &all_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (my_id == 0) {
        printf("%f\n", all_sum);
    }
}

int main(int argc, char *argv[]) {
    int my_id = 0;
    int size = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
//    MPI_Send(&sum, 1, MPI_INT, 0, 99, MPI_COMM_WORLD);
//    MPI_Recv(&receive, 1, MPI_INT, 1, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//    MPI_Finalize();

    const int n = 1000;
    calc_e(my_id, size, n);

    MPI_Finalize();
    return 0;
}
