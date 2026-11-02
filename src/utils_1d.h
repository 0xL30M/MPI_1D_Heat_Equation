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

#ifndef UTILS_1D_H
#define UTILS_1D_H

#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#ifdef USE_OMP
#include <omp.h>
#endif
#ifdef USE_MPI
#include <mpi.h>
#endif

#define ITER 20000

/**
  * @brief Parses the input file and extracts the positions and temperatures of each point.
  *
  * @param input_file Path to the input file containing the initial positions and temperatures.
  * @param x_pos Array of `double` values representing the x-coordinates of the bodies.
  * @param x_temp Array of `double` values representing the temperatures of the bodies.
  * @param num_pts Number of bodies in the simulation.
  *
  * @return 0 on success, -1 on error (e.g., file not found or inaccessible or wrongly formatted).
  */
static inline int parse_input(const char *input_file, double *x_pos, double *x_temp, size_t num_pts) {
  FILE* file = fopen(input_file, "r");
  if (!file) {
    fprintf(stderr, "Error opening file");
    return -1;
  }
  
  // Save the positions and temperatures of each body
  for (size_t i = 0; i < num_pts + 2; ++i) {
    if (fscanf(file, "%lf %lf", &x_pos[i], &x_temp[i]) != 2) {
      fprintf(stderr, "Error reading body %zu", i);
      fclose(file);
      return -1;
    }
  }

  return 0;
}

/**
 * @brief Writes the results of the simulation to an output file.
 *
 * @param output_file Path to the output file where the results will be written.
* @param err Error value of the simulation.
 * @param x Array of `double` values representing the x-coordinates of the points.
 * @param u Array of `double` values representing the temperatures of the points.
 * @param elapsed Time taken to complete the simulation.
 *
* @return 0 on success, -1 on error.
 */
static inline int write_solution(const char* output_file, double err, const double *x, const double *u,
                                  double elapsed, size_t num_pts)
{
  FILE *file;
  if ((file = fopen(output_file, "w")) == NULL) {
    fprintf(stderr, "Error opening file %s\n", output_file);
    return -1;
  }

  // Write the error
  fprintf(file, "%.4e\n", err);

  // Write the positions and temperatures
  for (size_t i = 0; i < num_pts + 2; ++i) {
    fprintf(file, "%.4e  %.4e\n", x[i], u[i]);
  }

  // Write the elapsed time
  fprintf(file, "%lf\n", elapsed);

  fclose (file);
  return 0;
}


int heat_eqn_solver_1d(double *heat_array, double *err, size_t num_pts);

#endif // UTILS_H
