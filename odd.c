#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(a,b)  ((a)<(b)?(a):(b))
#define c(x) printf("%s\n", #x)
#define e(x) printf("%s = %d\n", #x, x)

int main (int argc, char *argv[])
{
   int    count;        /* Local prime count */
   double elapsed_time; /* Parallel execution time */
   int    first;        /* Index of first multiple */
   int    global_count; /* Global prime count */
   int    high_value;   /* Highest value on this proc */
   int    i;
   int    id;           /* Process ID number */
   int    index;        /* Index of current prime */
   int    low_value;    /* Lowest value on this proc */
   char  *marked;       /* Portion of 2,...,'n' */
   int    n;            /* Sieving from 2, ..., 'n' */
   int    p;            /* Number of processes */
   int    proc0_size;   /* Size of proc 0's subarray */
   int    prime;        /* Current prime */
   int    size;         /* Elements in 'marked' */
   int    odd_size;

   MPI_Init (&argc, &argv);

   /* Start the timer */

   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();

   if (argc != 2) {
      if (!id) printf ("Command line: %s <m>\n", argv[0]);
      MPI_Finalize();
      exit (1);
   }

   n = atoi(argv[1]);

   /* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */

   low_value = 2 + id*(n-1)/p;
   if(low_value % 2 == 0) low_value++;
   high_value = 1 + (id+1)*(n-1)/p;
   size = high_value - low_value + 1;
   
   /* Bail out if all the primes used for sieving are
      not all held by process 0 */

   proc0_size = (n-1)/p;

   if ((2 + proc0_size) < (int) sqrt((double) n)) {
      if (!id) printf ("Too many processes\n");
      MPI_Finalize();
      exit (1);
   }

   if(size % 2 == 0) odd_size = size / 2;
   else odd_size = (size + 1) / 2;
   

   /* Allocate this process's share of the array. */

   marked = (char *) malloc (odd_size);

   if (marked == NULL) {
      printf ("Cannot allocate enough memory\n");
      MPI_Finalize();
      exit (1);
   }

   // memset
   for (i = 0; i < odd_size; i++) marked[i] = 0;
   
   if (!id) index = 0;

   prime = 3;
   do {
      // finding "first"
      if (prime * prime > low_value){
         first = (prime * prime - low_value) / 2;
      }
      else {
         if (!(low_value % prime)) first = 0;
         else{
            if((low_value % prime) % 2 != 0)
               first = (prime - (low_value % prime)) / 2;
            else first = prime - (low_value % prime) / 2;
         } 
      }

      //  sieving
      for (i = first; i < odd_size; i += prime) marked[i] = 1;
      
      if (!id) {
         while (marked[++index]);
         prime = index * 2 + 3;
      }
      
      if (p > 1) MPI_Bcast (&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
   } while (prime * prime <= n);

   // counting
   count = 0;
   for (i = 0; i < odd_size; i++)
      if (!marked[i]) count++;
   if (p > 1) MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM,
      0, MPI_COMM_WORLD);
   else global_count = count;

   /* Stop the timer */

   elapsed_time += MPI_Wtime();


   /* Print the results */

   if (!id) {
      // + 1 cz' 2 is not counted
      printf ("There are %d primes less than or equal to %d\n",
         global_count + 1, n);
      printf ("SIEVE (%d) %10.6f\n", p, elapsed_time);
   }
   MPI_Finalize ();
   return 0;
}