#include "utils_1d.h"

/* Implement the MPI solver here. You may add includes and helper functions.
 * Keep the explicit finite difference scheme and the public signature.
 * MPI is initialized/finalized by the driver. The full result and cumulative
 * error must be available on rank 0 when this function returns.
 */

 /**
 * @brief Solve the heat equation using an explicit finite difference scheme.
 *
 * This functions solves the heat equation iterating over ITER time steps.
 *
 * @param heat_array Array of `double` values representing the temperatures of the points.
 * @param err Pointer to a `double` value representing the error in the simulation.
 * @param num_pts Number of points in the simulation.
 *
 * @return 0 on success, -1 on error.
 */
int heat_eqn_solver_1d(double *heat_array, double *err, size_t num_pts){
  // Setup the parameters
  double dx = 1.0/(num_pts+1);
  double dt = 0.5*dx*dx;
  size_t beg = 1;
  size_t end = beg + num_pts - 1;

  // Utilize an auxiliary array to store the new temperatures
  double *aux = (double *) malloc((num_pts + 2) * sizeof(double));

  for (size_t i = 0; i < ITER; i++){

    // Set boundary conditions to 0
    aux[beg - 1] = 0.0;
    aux[end + 1] = 0.0;

    // Compute temperatures at new time step
    for (size_t j = beg; j <= end; j++){
      aux[j] = heat_array[j] + dt/(dx*dx)*(heat_array[j-1] - 2.0*heat_array[j] + heat_array[j+1]);
      *err += dx*fabs(aux[j] - heat_array[j]);
    }

    // Copy back into the original array
    for (size_t j = beg; j <= end; j++){
      heat_array[j] = aux[j];
    }
  }
  free(aux);
  return 0;
}
