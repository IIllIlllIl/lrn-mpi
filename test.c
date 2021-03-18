#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char* argv[]) {
    int done = 0, n, myid, numprocs, i;

    double mypi, pi, sum;
    double startwtime, endwtime;
    int  namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Get_processor_name(processor_name, &namelen);

    fprintf(stderr, "Process %d on %s\n", myid, processor_name);
    fflush(stderr);

    n = 0;
    while (!done) {
        if (myid == 0) {
            printf("输入一个数字不超过100000000: (0 退出) "); fflush(stdout);
            scanf("%d", &n);

            startwtime = MPI_Wtime();
        }
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);/*将n值广播出去*/
        if (n == 0)
            done = 1;
        else {

            sum = 0.0;
            for (i = myid + 1; i <= n; i += numprocs) {


                sum += i;
            }
            mypi = sum;/*各个进程并行计算得到的部分和*/

            MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            if (myid == 0) {
                /*执行累加的0号进程将近似值打印出来*/
                printf("结果 %.16f\n", pi);
                endwtime = MPI_Wtime();
                printf("时间 = %f\n", endwtime - startwtime);
            }
        }
    }
    MPI_Finalize();
}