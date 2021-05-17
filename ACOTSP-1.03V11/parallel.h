#include <mpi.h>

/* NUMERO DE PROCESOS TOTALES*/
int numprocs;
/* RANK */
int rank;
/* ESTADO DE UNA OPERACION DE RECEPCION*/
MPI_Status status;
/* LA MEJOR RUTA ENCONTRADA HASTA EL MOMENTO EN GLOBAL*/
long int best_global_tour;

void write_mpi_report ( void );
void BcastBestSolutionToColonies ( void );
void foreign_solution_update_pheromone( long int *ftour );
void foreign_solution_update_pheromone_weighted( long int *ftour, long int weight );



