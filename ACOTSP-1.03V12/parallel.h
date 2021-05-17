#include <mpi.h>

/* NUMERO DE PROCESOS TOTALES*/
int numprocs;
/* RANK */
int rank;
/* ESTADO DE UNA OPERACION DE RECEPCION*/
MPI_Status status;

void write_mpi_report ( void );
/*void BcastBestSolutionToColonies ( void );*/
void startCommColonies ( void );
void sendBestSolutionToColonies ( void );
void listenColonies( void );

void foreign_solution_update_pheromone( long int *ftour );
void foreign_solution_update_pheromone_weighted( long int *ftour, long int weight );

