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
    int rank, size;
    
    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) return -1;
    if (MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS) return -1;

    size_t rem = num_pts % size;
    size_t local_n = (num_pts / size) + ((size_t)rank < rem ? 1 : 0);
    size_t offset = ((size_t)rank < rem) ? (rank * local_n) : (rank * (num_pts / size) + rem);

    double dx = 1.0 / (num_pts + 1);

    // Allocazione allineata per l'uso dei registri AVX-512
    double *buffer;
    size_t alloc_size = 2 * (local_n + 2) * sizeof(double);
    if (posix_memalign((void**)&buffer, 64, alloc_size) != 0) return -1;

    // Inizializzazione a 0
    for(size_t i = 0; i < 2 * (local_n + 2); i++) {
        buffer[i] = 0.0;
    }

    double * restrict curr = buffer;
    double * restrict next = buffer + (local_n + 2);

    if (local_n > 0) {
        for (size_t i = 0; i < local_n; i++) {
            curr[i + 1] = heat_array[offset + 1 + i];
        }
    }

    // Ricerca dei vicini (Topological search anti-deadlock)
    int left_neighbor = MPI_PROC_NULL;
    int right_neighbor = MPI_PROC_NULL;

    if (local_n > 0) {
        for (int r = rank - 1; r >= 0; r--) {
            if ((num_pts / size) + ((size_t)r < rem ? 1 : 0) > 0) { left_neighbor = r; break; }
        }
        for (int r = rank + 1; r < size; r++) {
            if ((num_pts / size) + ((size_t)r < rem ? 1 : 0) > 0) { right_neighbor = r; break; }
        }
    }

    double step_err = 0.0;

    // Loop Temporale Principale
    for (size_t i = 0; i < ITER; i++){
        if (local_n > 0) {
            // Halo Exchange bloccante contiguo (perfetto per latenze bassissime < 1 byte)
            MPI_Sendrecv(&curr[1], 1, MPI_DOUBLE, left_neighbor, 0,
                         &curr[local_n + 1], 1, MPI_DOUBLE, right_neighbor, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Sendrecv(&curr[local_n], 1, MPI_DOUBLE, right_neighbor, 1,
                         &curr[0], 1, MPI_DOUBLE, left_neighbor, 1,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            double iter_err = 0.0;
            
            // Loop computazionale unico e linearizzato.
            // Grazie alla matematica semplificata, richiede la metà delle operazioni CPU.
            #pragma GCC ivdep
            for (size_t j = 1; j <= local_n; j++){
                double new_val = 0.5 * (curr[j-1] + curr[j+1]);
                iter_err += fabs(new_val - curr[j]);
                next[j] = new_val;
            }

            step_err += iter_err;
        }
        
        // Pointer Swap
        double * restrict temp = curr;
        curr = next;
        next = temp;
    }

    int *recvcounts = NULL;
    int *displs = NULL;

    if (rank == 0) {
        recvcounts = (int *)malloc(size * sizeof(int));
        displs = (int *)malloc(size * sizeof(int));
        
        if (!recvcounts || !displs) {
            if(recvcounts) free(recvcounts);
            if(displs) free(displs);
            free(buffer);
            return -1;
        }
        
        for (int r = 0; r < size; r++) {
            size_t r_n = (num_pts / size) + ((size_t)r < rem ? 1 : 0);
            size_t r_off = ((size_t)r < rem) ? (r * r_n) : (r * (num_pts / size) + rem);
            recvcounts[r] = (int)r_n;
            displs[r] = (int)(r_off + 1);
        }
    }

    MPI_Gatherv(&curr[1], (int)local_n, MPI_DOUBLE,
                heat_array, recvcounts, displs, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    // Riduzione finale dell'errore (moltiplicato per dx solo una volta)
    double local_err = step_err * dx; 
    double global_err = 0.0;
    MPI_Reduce(&local_err, &global_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        *err += global_err;
        heat_array[0] = 0.0;
        heat_array[num_pts + 1] = 0.0;
        
        free(recvcounts);
        free(displs);
    }

    free(buffer);
    return 0;
}
