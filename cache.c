#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define cache 12 * 1024 * 1024

int main(int argc, char *argv[])
{
   int count;           /* Local prime count */
   double elapsed_time; /* Parallel execution time */
   int first;           /* Index of first multiple */
   int global_count;    /* Global prime count */
   int i;
   int j;
   int id;         /* Process ID number */
   int index;      /* Index of current prime */
   int low_value;  /* Lowest value on this proc */
   int high_value; /* Highest value on this proc */
   int low;
   int high;
   char *marked;   /* Portion of 2,...,'n' */
   int n;          /* Sieving from 2, ..., 'n' */
   int p;          /* Number of processes */
   int proc0_size; /* Size of proc 0's subarray */
   int prime;      /* Current prime */
   int *prm;
   char *visit;
   int size; /* Elements in 'marked' */
   int odd_size;
   int block;

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
   low_value = 2 + id * ((n - 1) / p);
   if (low_value % 2 == 0)
      low_value++;
   high_value = 1 + (id + 1) * ((n - 1) / p);
   if (high_value % 2 != 0)
      high_value++;
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
   block = (cache / sizeof(int)) / p;
   marked = (char *)malloc(block);
   if (marked == NULL)
   {
      printf("Cannot allocate enough memory\n");
      MPI_Finalize();
      exit(1);
   }

   // prime in 0 ~ sqrt(n)
   visit = (char *)malloc((int)sqrt((double)n));
   prm = (int *)malloc((int)sqrt((double)n) * sizeof(int));

   for (i = 0; i < (int)sqrt((double)n); i++)
      visit[i] = 0;
   for (i = 0; i < (int)sqrt((double)n); i++)
      prm[i] = 0;

   for (i = 2; i <= (int)sqrt((double)n); i++)
   {
      if (!visit[i])
      {
         prm[++prm[0]] = i;
      }
      for (j = 1; j <= prm[0] && i * prm[j] <= (int)sqrt((double)n); j++)
      {
         visit[i * prm[j]] = 1;
         if (i % prm[j] == 0)
         {
            break;
         }
      }
   }

   count = 0;

   for (low = low_value; low < high_value; low += (block << 1))
   {
      high = MIN(low + (block << 1) - 1, high_value);

      // memset
      for (i = 0; i < block; i++)
      {
         marked[i] = 0;
      }

      index = 2;
      do
      {
         prime = prm[index++];
         if (prime == 0)
            break;

         // finding "first"
         if (prime > (int)sqrt((double)low))
         {
            first = (prime * prime - low) >> 1;
         }
         else
         {
            if (!(low % prime))
               first = 0;
            else
            {
               if ((low % prime) % 2 != 0)
                  first = (prime - (low % prime)) >> 1;
               else
                  first = prime - ((low % prime) >> 1);
            }
         }

         //  sieving
         for (i = first; i < ((high - low + 1) >> 1); i += prime)
            marked[i] = 1;
      } while (prime <= (int)sqrt((double)high));

      // counting
      for (i = 0; i < ((high - low + 1) >> 1); i++)
         if (!marked[i])
         {
            count++;
         }
   }

   if (p > 1)
      MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM,
                 0, MPI_COMM_WORLD);
   else
      global_count = count;

   /* Stop the timer */

   elapsed_time += MPI_Wtime();

   /* Print the results */

   if (!id)
   {
      // + 1 cz' 2 is not counted
      printf("There are %d primes less than or equal to %d\n",
             global_count + 1, n);
      printf("SIEVE (%d) %10.6f\n", p, elapsed_time);
   }
   MPI_Finalize();
   return 0;
}