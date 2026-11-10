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

    // Distribuzione bilanciata del carico di lavoro
    size_t rem = num_pts % size;
    size_t local_n = (num_pts / size) + ((size_t)rank < rem ? 1 : 0);
    size_t offset = ((size_t)rank < rem) ? (rank * local_n) : (rank * (num_pts / size) + rem);

    double dx = 1.0 / (num_pts + 1);
    double dt = 0.5 * dx * dx;
    const double factor = dt / (dx * dx); // Pre-calcolato per costare zero nei loop

    // Allineamento a 64 bytes per massimizzare il throughput delle istruzioni AVX-512
    double *buffer;
    size_t alloc_size = 2 * (local_n + 2) * sizeof(double);
    if (posix_memalign((void**)&buffer, 64, alloc_size) != 0) return -1;

    // Inizializzazione manuale (evita garbage nelle celle fantasma ai bordi esterni)
    for(size_t i = 0; i < 2 * (local_n + 2); i++) {
        buffer[i] = 0.0;
    }

    // restrict: informa il compilatore che i puntatori non aliasano memoria condivisa
    double * restrict curr = buffer;
    double * restrict next = buffer + (local_n + 2);

    if (local_n > 0) {
        for (size_t i = 0; i < local_n; i++) {
            curr[i + 1] = heat_array[offset + 1 + i];
        }
    }

    // Topological Neighbor Search (Previene Deadlock in scenari con più Rank che Punti)
    int left_neighbor = MPI_PROC_NULL;
    int right_neighbor = MPI_PROC_NULL;

    if (local_n > 0) {
        for (int r = rank - 1; r >= 0; r--) {
            size_t r_n = (num_pts / size) + ((size_t)r < rem ? 1 : 0);
            if (r_n > 0) { left_neighbor = r; break; }
        }
        for (int r = rank + 1; r < size; r++) {
            size_t r_n = (num_pts / size) + ((size_t)r < rem ? 1 : 0);
            if (r_n > 0) { right_neighbor = r; break; }
        }
    }

    double step_err = 0.0;

    // Loop Temporale Principale
    for (size_t i = 0; i < ITER; i++){
        if (local_n > 0) {
            MPI_Request reqs[4];
            int req_cnt = 0;
            
            // 1. COMUNICAZIONE ASINCRONA: Inizio trasferimento Halo/Ghost Cells
            if (left_neighbor != MPI_PROC_NULL) {
                MPI_Irecv(&curr[0], 1, MPI_DOUBLE, left_neighbor, 1, MPI_COMM_WORLD, &reqs[req_cnt++]);
                MPI_Isend(&curr[1], 1, MPI_DOUBLE, left_neighbor, 0, MPI_COMM_WORLD, &reqs[req_cnt++]);
            }
            if (right_neighbor != MPI_PROC_NULL) {
                MPI_Irecv(&curr[local_n + 1], 1, MPI_DOUBLE, right_neighbor, 0, MPI_COMM_WORLD, &reqs[req_cnt++]);
                MPI_Isend(&curr[local_n], 1, MPI_DOUBLE, right_neighbor, 1, MPI_COMM_WORLD, &reqs[req_cnt++]);
            }

            double iter_err = 0.0;
            
            // 2. COMPUTAZIONE SOVRAPPOSTA: Calcoliamo il dominio INTERNO.
            // Mentre l'hardware di rete lavora, la CPU macina dati in parallelo.
            #pragma GCC ivdep
            for (size_t j = 2; j < local_n; j++){
                double diff = factor * (curr[j-1] - 2.0*curr[j] + curr[j+1]);
                next[j] = curr[j] + diff;
                iter_err += fabs(diff);
            }

            // 3. SINCRONIZZAZIONE: Ora dobbiamo aspettare rigorosamente i bordi
            if (req_cnt > 0) {
                MPI_Waitall(req_cnt, reqs, MPI_STATUSES_IGNORE);
            }

            // 4. COMPUTAZIONE DEI BORDI: calcolati in sicurezza dopo la Wait
            if (local_n == 1) {
                double diff = factor * (curr[0] - 2.0*curr[1] + curr[2]);
                next[1] = curr[1] + diff;
                iter_err += fabs(diff);
            } else if (local_n > 1) {
                // Bordo Sinistro
                double diff_l = factor * (curr[0] - 2.0*curr[1] + curr[2]);
                next[1] = curr[1] + diff_l;
                iter_err += fabs(diff_l);
                
                // Bordo Destro
                double diff_r = factor * (curr[local_n-1] - 2.0*curr[local_n] + curr[local_n+1]);
                next[local_n] = curr[local_n] + diff_r;
                iter_err += fabs(diff_r);
            }

            step_err += iter_err;
        }
        
        // Pointer Swap a Costo Zero (invece di ricopiare l'array byte per byte)
        double * restrict temp = curr;
        curr = next;
        next = temp;
    }

    // Set-up dinamico per Gatherv
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

    // Ricostruzione dell'array master sul Rank 0
    MPI_Gatherv(&curr[1], (int)local_n, MPI_DOUBLE,
                heat_array, recvcounts, displs, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    // Ottimizzazione matematica: L'errore viene moltiplicato per dx solo una volta alla fine
    double local_err = step_err * dx; 
    double global_err = 0.0;
    
    // Riduzione finale degli errori
    MPI_Reduce(&local_err, &global_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        *err += global_err;
        // Reinserimento forzato dei bordi statici
        heat_array[0] = 0.0;
        heat_array[num_pts + 1] = 0.0;
        
        free(recvcounts);
        free(displs);
    }

    free(buffer);
    return 0;
}
