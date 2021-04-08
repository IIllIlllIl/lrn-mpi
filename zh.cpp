#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <iostream>
#define MIN(a, b) ((a) < (b) ? (a) : (b))

int main(int argc, char *argv[])
{
    int count;           /* Local prime count */
    double elapsed_time; /* Parallel execution time */
    int first;           /* Index of first multiple */
    int global_count;    /* Global prime count */
    int high_value;      /* Highest value on this proc */
    int i;
    int id;        /* Process ID number */
    int index;     /* Index of current prime */
    int low_value; /* Lowest value on this proc */
    char *marked;  /* Portion of 2,...,'n' */
    char *pmarked;
    int n;          /* Sieving from 2, ..., 'n' */
    int p;          /* Number of processes */
    int proc0_size; /* Size of proc 0's subarray */
    int prime;      /* Current prime */
    int size;       /* Elements in 'marked' */
    int flag_over = 0;

    MPI_Init(&argc, &argv);

    /* Start the timer */

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    if (argc != 2)
    {
        if (!id)
            printf("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    n = atoi(argv[1]);

    /* Figure out this process's share of the array, as
       well as the integers represented by the first and
       last array elements */

    low_value = 2 + id * (n - 1) / p;
    high_value = 1 + (id + 1) * (n - 1) / p;
    size = high_value - low_value + 1;

    /* Bail out if all the primes used for sieving are
       not all held by process 0 */

    proc0_size = (n - 1) / p;

    if ((2 + proc0_size) < (int)sqrt((double)n))
    {
        if (!id)
            printf("Too many processes\n");
        MPI_Finalize();
        exit(1);
    }
    /* Allocate this process's share of the array. */
    int sqn = sqrt(n);

    marked = (char *)malloc(size);
    pmarked = (char *)malloc(size);
    for (i = 0; i < size; i++)
        marked[i] = pmarked[i] = 0;
    bool setflag = low_value % 2;

    prime = 3;
    index = 1;
    do
    {
        for (int i = prime * prime - 2; i < size; i += prime)
            if (i % 2 ^ setflag)
                pmarked[i] = 1;
        while (pmarked[((++index) << 1) - 1])
            ;
        prime = (index << 1) + 1;
    } while (prime * prime <= sqn);

    index = 1;
    prime = 3;
    do
    {
        if (prime * prime > low_value)
            first = prime * prime - low_value;
        else
        {
            if (!(low_value % prime))
                first = 0;
            else
                first = prime - (low_value % prime);
        }
        for (int i = first; i < size; i += prime)
            if (i % 2 ^ setflag)
                marked[i] = 1;
        while (pmarked[((++index) << 1) - 1])
            ;
        prime = (index << 1) + 1;
    } while (prime <= sqn);

    count = 0;
    for (i = !(low_value % 2); i < size; i += 2)
        if (!marked[i])
            count++;
    if (p > 1)
        MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM,
                   0, MPI_COMM_WORLD);

    /* Stop the timer */

    elapsed_time += MPI_Wtime();

    /* Print the results */

    if (!id)
    {
        printf("There are %d primes less than or equal to %d\n",
               global_count + 1, n);
        printf("SIEVE (%d) %10.6f\n", p, elapsed_time);
    }
    MPI_Finalize();
    return 0;
}