#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// size of cache
#define cache 12 * 1024 * 1024

int main(int argc, char *argv[])
{
   int count;           /* Local prime count */
   double elapsed_time; /* Parallel execution time */
   int first;           /* Index of first multiple */
   int global_count;    /* Global prime count */
   int i;
   int j;
   int id;   /* Process ID number */
   int *ind; /* Index of current prime */
   int index;
   int low_value;  /* Lowest value on this proc */
   int high_value; /* Highest value on this proc */
   int low;
   int high;
   char *marked;   /* Portion of 2,...,'n' */
   int n;          /* Sieving from 2, ..., 'n' */
   int p;          /* Number of processes */
   int proc0_size; /* Size of proc 0's subarray */
   int prime;      /* Current prime */
   int prime_mult_2;
   int *prm;
   char *visit;
   int size; /* Elements in 'marked' */
   int odd_size;
   int block;
   int sqrt_n;
   int sqrt_high;
   int block_size;

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
   sqrt_n = (int)sqrt((double)n);
   if (n > 1000)
   {
      /* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */
      low_value = 2 + (long long)id * (n - 1) / p;
      if (!(low_value & 0x1))
         low_value++;
      high_value = 1 + (long long)(id + 1) * (n - 1) / p;
      if (high_value & 0x1)
         high_value++;
      size = high_value - low_value + 1;

      /* Bail out if all the primes used for sieving are
      not all held by process 0 */
      proc0_size = (n - 1) / p;
      if ((2 + proc0_size) < sqrt_n)
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

      sqrt_n++;

      // prime in 0 ~ sqrt(n)
      visit = (char *)malloc(sqrt_n);
      prm = (int *)malloc(sqrt_n * sizeof(int));

      memset(visit, 0, sqrt_n);
      memset(prm, 0, sqrt_n * sizeof(int));

      for (i = 2; i <= sqrt_n; i++)
      {
         if (!visit[i])
         {
            prm[++prm[0]] = i;
         }
         for (j = 1; j <= prm[0] && i * prm[j] <= sqrt_n; j++)
         {
            visit[i * prm[j]] = 1;
            if (!(i % prm[j]))
            {
               break;
            }
         }
      }

      count = 0;

      for (low = low_value; low < high_value; low += (block << 1))
      {
         high = (low + (block << 1) - 1) < high_value ? (low + (block << 1) - 1) : high_value;
         sqrt_high = (int)sqrt((double)high);
         block_size = (high - low + 1) >> 1;

         // memset
         memset(marked, 0, block);

         ind = prm + 2;
         do
         {
            prime = *ind++;
            if (!prime)
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
                  if ((low % prime) & 0x1)
                     first = (prime - (low % prime)) >> 1;
                  else
                     first = prime - ((low % prime) >> 1);
               }
            }

            prime_mult_2 = prime << 1;
            //  sieving
            for (i = first; i < block_size; i += prime_mult_2 << 1)
            {
               marked[i] = 1;
               if (i + prime < block_size)
               {
                  marked[i + prime] = 1;
               }
               if (i + prime_mult_2 < block_size)
               {
                  marked[i + prime_mult_2] = 1;
               }
               if (i + prime + prime_mult_2 < block_size)
               {
                  marked[i + prime + prime_mult_2] = 1;
               }
            }
         } while (prime <= sqrt_high);

         // counting
         for (i = 0; i < (block_size) >> 1; i += 2)
         {
            if (!marked[i])
            {
               count++;
            }
            if (!marked[i + 1])
            {
               count++;
            }
            if (!marked[block_size - i - 1])
            {
               count++;
            }
            if (block_size - (i << 1) - 3)
               if (!marked[block_size - i - 2])
               {
                  count++;
               }
         }
      }
   }
   else
   {
      /* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */

      low_value = 2 + id * (n - 1) / p;
      if (!(low_value & 0x1))
         low_value++;
      high_value = 1 + (id + 1) * (n - 1) / p;
      size = high_value - low_value + 1;

      /* Bail out if all the primes used for sieving are
      not all held by process 0 */

      proc0_size = (n - 1) / p;

      if ((2 + proc0_size) < sqrt_n)
      {
         if (!id)
            printf("Too many processes\n");
         MPI_Finalize();
         exit(1);
      }

      if (!(size & 0x1))
         odd_size = size >> 1;
      else
         odd_size = (size + 1) >> 1;

      /* Allocate this process's share of the array. */

      marked = (char *)malloc(odd_size);

      if (marked == NULL)
      {
         printf("Cannot allocate enough memory\n");
         MPI_Finalize();
         exit(1);
      }

      // memset
      for (i = 0; i < odd_size; i++)
         marked[i] = 0;

      if (!id)
         index = 0;

      prime = 3;
      do
      {
         // finding "first"
         if (prime * prime > low_value)
         {
            first = (prime * prime - low_value) >> 1;
         }
         else
         {
            if (!(low_value % prime))
               first = 0;
            else
            {
               if ((low_value % prime) & 0x1)
                  first = (prime - (low_value % prime)) >> 1;
               else
                  first = prime - ((low_value % prime) >> 1);
            }
         }

         //  sieving
         prime_mult_2 = prime << 1;
         for (i = first; i < odd_size; i += prime_mult_2 << 1)
         {
            marked[i] = 1;
            if (i + prime < odd_size)
            {
               marked[i + prime] = 1;
            }
            if (i + prime_mult_2 < odd_size)
            {
               marked[i + prime_mult_2] = 1;
            }
            if (i + prime + prime_mult_2 < odd_size)
            {
               marked[i + prime + prime_mult_2] = 1;
            }
         }

         if (!id)
         {
            while (marked[++index])
               ;
            prime = (index << 1) + 3;
         }

         if (p > 1)
            MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
      } while (prime * prime <= n);

      // counting
      count = 0;
      for (i = 0; i < (odd_size) >> 1; i += 2)
      {
         if (!marked[i])
         {
            count++;
         }
         if (!marked[i + 1])
         {
            count++;
         }
         if (!marked[odd_size - i - 1])
         {
            count++;
         }
         if (odd_size - (i << 1) - 3)
            if (!marked[odd_size - i - 2])
            {
               count++;
            }
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
