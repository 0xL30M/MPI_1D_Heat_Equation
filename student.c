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

    // Padding estremo: allineamento forzato al confine della cache line (8 double = 64 byte).
    // Questo previene i Cache Line Split garantendo che c[8] e n[8] siano 64-byte aligned.
    size_t stride = (local_n + 15) & ~7; // Arrotondamento al multiplo di 8 superiore
    size_t alloc_size = 2 * (stride + 16) * sizeof(double);
    
    double *buffer;
    if (posix_memalign((void**)&buffer, 64, alloc_size) != 0) return -1;

    for(size_t i = 0; i < 2 * (stride + 16); i++) {
        buffer[i] = 0.0;
    }

    double * restrict curr = buffer;
    double * restrict next = buffer + stride + 16;

    // offset + 8 per saltare la regione fantasma allargata usata per l'allineamento
    if (local_n > 0) {
        for (size_t i = 0; i < local_n; i++) {
            curr[8 + i] = heat_array[offset + 1 + i];
        }
    }

    // Topological search anti-deadlock
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

    // Inizializzazione MPI Persistent Communications con Double Buffering
    MPI_Request reqs[2][4];
    int req_counts[2] = {0, 0};

    for (int step = 0; step < 2; step++) {
        double *c = (step == 0) ? curr : next;
        
        if (left_neighbor != MPI_PROC_NULL) {
            MPI_Recv_init(&c[7], 1, MPI_DOUBLE, left_neighbor, 1, MPI_COMM_WORLD, &reqs[step][req_counts[step]++]);
            MPI_Send_init(&c[8], 1, MPI_DOUBLE, left_neighbor, 0, MPI_COMM_WORLD, &reqs[step][req_counts[step]++]);
        }
        if (right_neighbor != MPI_PROC_NULL) {
            MPI_Recv_init(&c[8 + local_n], 1, MPI_DOUBLE, right_neighbor, 0, MPI_COMM_WORLD, &reqs[step][req_counts[step]++]);
            MPI_Send_init(&c[8 + local_n - 1], 1, MPI_DOUBLE, right_neighbor, 1, MPI_COMM_WORLD, &reqs[step][req_counts[step]++]);
        }
    }

    double step_err = 0.0;

    // Loop Temporale Principale
    for (size_t i = 0; i < ITER; i++){
        if (local_n > 0) {
            int step_idx = i & 1; // Equivalente a i % 2, ottimizzato a bitwise
            
            // Startall elimina la latenza di attraversamento dello stack rispetto a Sendrecv
            if (req_counts[step_idx] > 0) {
                MPI_Startall(req_counts[step_idx], reqs[step_idx]);
                MPI_Waitall(req_counts[step_idx], reqs[step_idx], MPI_STATUSES_IGNORE);
            }

            double iter_err = 0.0;
            
            // Assunzione esplicita di allineamento per vettorializzazione AVX-512 perfetta senza peeled loops
            double * restrict c = __builtin_assume_aligned(curr, 64);
            double * restrict n = __builtin_assume_aligned(next, 64);
            
            #pragma GCC ivdep
            for (size_t j = 8; j < local_n + 8; j++){
                double new_val = 0.5 * (c[j-1] + c[j+1]);
                iter_err += fabs(new_val - c[j]);
                n[j] = new_val;
            }

            step_err += iter_err;
        }
        
        // Zero-copy pointer swap
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

    MPI_Gatherv(&curr[8], (int)local_n, MPI_DOUBLE,
                heat_array, recvcounts, displs, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    // Riduzione finale globale per l'errore calcolato
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

    // Disallocazione risorse MPI persistenti
    for (int step = 0; step < 2; step++) {
        for (int k = 0; k < req_counts[step]; k++) {
            MPI_Request_free(&reqs[step][k]);
        }
    }

    free(buffer);
    return 0;
}
