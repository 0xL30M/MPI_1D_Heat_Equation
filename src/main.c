/**
 * **************************************************************
 * WARNING: THIS FILE IS NOT MODIFIABLE UNDER ANY CIRCUMSTANCES.
 * **************************************************************
 *
 * Any modifications to this file will result in the 
 * INVALIDATION of your submission. During testing, we will verify 
 * that this file has not been altered in any way.
 * 
 * Please adhere strictly to the original structure and code in 
 * this file. Any changes made to this file will be detected, and 
 * the run will be considered invalid.
 *
 * **************************************************************
 */

#include "utils_1d.h"

/**
 * Setup macros based on the compilation flags.
 * The macro set up are:
 * __RANK_0__ : Execute the expression only if the rank is 0 or if not using MPI
 * __TIME_START__ : Start the timer
 * __TIME_END__ : End the timer
 * __TIME_MSECS__ : Calculate the elapsed time in milliseconds, for MPI do a reduction
 *
 * Also start time, end time and elapsed time variables are declared
 * with different types based on the compilation flags
 **/

#undef USE_OMP

#ifndef USE_MPI
  #define USE_MPI
#endif

#if defined(USE_MPI)
#define __RANK_0__(expr)  do{if (rank == 0){expr}}while(0)
#define __TIME_START__    (time_start    = MPI_Wtime())
#define __TIME_END__      (time_end      = MPI_Wtime())
#define __TIME_MSECS__    do{                     \
                            elapsed = (time_end - time_start) * 1000.0; \
                            if (rank == 0) { \
                              ret = MPI_Reduce(MPI_IN_PLACE, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); \
                            } else { ret = MPI_Reduce(&elapsed, NULL, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); }\
                              if (MPI_SUCCESS != ret) goto err_hndl; \
                            }while(0)
double time_start, time_end, elapsed;
#endif

#if defined(USE_OMP) && !defined(USE_MPI)
#define __RANK_0__(expr)  do{expr}while(0)
#define __TIME_START__    (time_start    = omp_get_wtime())
#define __TIME_END__      (time_end      = omp_get_wtime())
#define __TIME_MSECS__    (elapsed = (time_end - time_start) * 1000.0)
double time_start, time_end, elapsed;
#endif

#if !defined(USE_OMP) && !defined(USE_MPI)
#include <time.h>
#define __RANK_0__(expr)  do{expr}while(0)
#define __TIME_START__    (time_start    = clock())
#define __TIME_END__      (time_end      = clock())
#define __TIME_MSECS__    (elapsed = (double)(((time_end - time_start) * 1000.0) / CLOCKS_PER_SEC ))
clock_t time_start, time_end;
double elapsed;
#endif // not defined USE_OMP && USE_MPI


int main(int argc, char ** argv)
{
  if (argc!=4){
    fprintf(stderr,"Usage: %s <num_pts> <input_file> <output_file>\n",argv[0]);
    return EXIT_FAILURE;
  }
  // Get number of points
  size_t num_pts = (size_t) atoll(argv[1]);
  if (num_pts <= 0){
    fprintf (stderr, "num_pts must be >0\n");
    return EXIT_FAILURE;
  }

  #ifdef USE_MPI
  MPI_Init(NULL,NULL);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  #endif

  
  int line = 0, ret;

  // Allocate memory for positions and temperatures
  double err = 0.0;
  double *x_pos = NULL, *heat_array = NULL;
  x_pos = (double *) malloc((num_pts + 2) * sizeof(double));
  heat_array = (double *) malloc((num_pts + 2) * sizeof(double));
  if (x_pos == NULL || heat_array == NULL){
    line = __LINE__;
    fprintf(stderr, "Error allocating x_pos or heat_array");
    goto err_hndl;
  }

  // Parse initial positions and temperatures from input file
  ret = parse_input(argv[2], x_pos, heat_array, num_pts);
  if (ret != 0) {line = __LINE__; goto err_hndl;}


  // Solve the heat equation
  __TIME_START__;
  ret = heat_eqn_solver_1d(heat_array, &err, num_pts);
  if (ret != 0) {line = __LINE__; goto err_hndl;}
  __TIME_END__;
  
  __TIME_MSECS__;

  // Write the solution to the output file
  __RANK_0__({
    ret = write_solution(argv[3], err, x_pos, heat_array, elapsed, num_pts);
    if (ret != 0) {line = __LINE__; goto err_hndl;}
    printf("Execution completed successfully in %f milliseconds.\n", elapsed);
  });

  // Free memory and finalize MPI if needed
  free(x_pos);
  free(heat_array);
#ifdef USE_MPI
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;

err_hndl: // Error handling
  fprintf(stderr, "Error at line %d\n", line);
  (void) line; // Silence compiler warning
  if (x_pos != NULL) free(x_pos);
  if (heat_array != NULL) free(heat_array);
#ifdef USE_MPI
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  return EXIT_FAILURE;
}

