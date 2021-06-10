/**
 
 EN UN FICHERO parallel.h PUEDES METER LAS VARIABLES GLOBALES QUE NECESITES PASAR AL 
 RESTO DE FUNCIONES. POR EJEMPLO: NPROC=NUMERO DE PROCESOS TOTALES, O ID=RANK DE CADA
 PROCESO, Y OTRAS COMO best_global_tour CON LA MEJOR RUTA ENCONTRADA HASTA EL MOMENTO
 EN GLOBAL...
 
 **/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <time.h>

#include "InOut.h"
#include "TSP.h"
#include "ants.h"
#include "ls.h"
#include "utilities.h"
#include "timer.h"
#include "parallel.h"


/***************************  Final MPI-Report ************************************/

/**
 Routine to write one summarize mpi report
 
 EMPEZAR POR ESTA RUTINA, PENSAR EN UN PROGRAMA MPI EMBARAZOSAMENTE PARALELO, TODOS
 LOS PROCESOS VAN POR LIBRE, PERO EN LA FUNCIÓN "exit_try" TODOS COMUNICAN LOS RESULTADOS 
 AL PROCESO 0 PARA QUE ESTE LOS IMPRIMA EN UN ÚNICO FICHERO, SIMILAR AL FICHERO "best.xxx"
 
 HAY QUE COMUNICAR: LONGITUD DEL MEJOR CAMINO, NÚMERO DE ITERACIONES, B-FAC, Y TIEMPO DE
 ENCONTRAR LA MEJOR SOLUCIÓN. 
 
 EL PROCESO 0 COMPRUEBA CUAL ES LA MEJOR SOLUCIÓN DE LAS QUE LE VAN LLEGANDO Y AL FINAL
 IMPRIME SOLO LA MEJOR DE ESE  "try"
 
 ACABA CON UN MPI_BARRIER() PARA QUE TODOS COMIENCEN A LA VEZ EL SIGUIENTE "try"
 **/
void write_mpi_report ( void )
{
    int rankSource, i;
    long tour_l, best_tour_length, iterations, best_iterations;
    double b_fact, best_b_fact, time, best_time;
        
    rankSource = rank;
    best_tour_length = best_so_far_ant->tour_length;
    best_iterations = found_best;
    best_b_fact = found_branching;
    best_time = time_used;
    
    if(rank == 0){

        for(i=1;i<numprocs;i++){
         MPI_Recv(&tour_l,1,MPI_LONG,i,1,MPI_COMM_WORLD,&status);
         MPI_Recv(&iterations,1,MPI_LONG,i,2,MPI_COMM_WORLD,&status);
         MPI_Recv(&b_fact,1,MPI_DOUBLE,i,3,MPI_COMM_WORLD,&status);
         MPI_Recv(&time,1,MPI_DOUBLE,i,4,MPI_COMM_WORLD,&status);
         
         if(tour_l < best_tour_length){
            rankSource = i;
            best_tour_length = tour_l;
            best_iterations = iterations;
            best_b_fact = b_fact;
            best_time = time;
         }
         
        }
        
        write_totalResults(rankSource, best_tour_length, best_iterations, best_b_fact,
              best_time);
       
    } else{
         MPI_Send(&best_tour_length,1,MPI_LONG,0,1,MPI_COMM_WORLD);
         MPI_Send(&best_iterations,1,MPI_LONG,0,2,MPI_COMM_WORLD);
         MPI_Send(&best_b_fact,1,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
         MPI_Send(&best_time,1,MPI_DOUBLE,0,4,MPI_COMM_WORLD);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    imprimir_resultados_parciales();

}



/******************************************************************/
/************************** Syncronous ****************************/

/*
 Bcast the best solution to the rest of the colonies.
 Called when a syncronous version is invoked
 
 UNA VEZ QUE LA VERSIÓN EMBARAZOSAMENTE PARALELA TE FUNCIONE. IMPLEMENTAMOS UNA SÍNCRONA.
 EN LA FUNCIÓN "control_xxx"  PODEMOS INCLUIR UN BCAST DE TODOS A TODOS COMUNICANDO SU MEJOR
 SOLUCIÓN HASTA ESE MOMENTO. SOLO VAN A COMUNICAR LA RUTA, NO LA HORMIGA COMPLETA, PORQUE 
 ASÍ SOLO HAY QUE COMUNICAR UN VECTOR DE LONG INTS DE N+1 ELEMENTOS. COMUNICAR LA HORMIGA
 HARÍA QUE TENGAMOS QUE COMUNICAR VARIOS MENSAJES O BIEN HACER UN NUEVO MPI_DATATYPE QUE ES 
 MÁS COMPLICADO, CUANDO LA ÚNICA INFORMACIÓN QUE QUEREMOS ES LA RUTA BUENA.
 
 UNA VEZ QUE RECIBEN LA MEJOR SOLUCIÓN DE CADA COLONIA TIENEN QUE LLAMAR A LA FUNCIÓN QUE 
 ACTUALIZA LAS FEROMONAS. ESA FUNCIÓN "update_pheromone" YA IMPLEMENTADA TIENE COMO ARGUMENTO
 LA HORMIGA. EN ESTE CASO COMO NO NOS PASAN LA HORMIGA SINO LA RUTA, PODEMOS HACER NUESTRA 
 PROPIA FUNCION DE ACTUALIZACIÓN DE FEROMONAS CON LAS SOLUCIONES EXTERNAS A LA QUE LE PASEMOS
 LA RUTA EN LUGAR DE LA HORMIGA COMPLETA.
 */
void BcastBestSolutionToColonies ( void )
{
  long int *tour_best_ant;
  long int distanceOwn, distanceForeing;
  int i, numBytes;
  
  tour_best_ant = calloc(n+1, sizeof(long int));
  distanceOwn =  compute_tour_length(best_so_far_ant->tour);
    
  for(i=0; i<numprocs; i++){
  
     if (i == rank){
     	numBytes= 8*(n+1);
        memcpy(tour_best_ant, best_so_far_ant->tour, numBytes);
       /* write_registro_antes(tour_best_ant, distanceOwn, best_so_far_ant->tour_length);*/
     }
         
     MPI_Bcast(tour_best_ant, n+1, MPI_LONG, i, MPI_COMM_WORLD);
     
     if (rank != i ){
        distanceForeing = compute_tour_length(tour_best_ant);
        /*write_registro_despues(tour_best_ant, distanceForeing);*/
        
        foreign_solution_update_pheromone(tour_best_ant);
        
        /*añadir lógica mejor camino que el propio*/
        if ( eas_flag && distanceOwn > distanceForeing){
           foreign_solution_update_pheromone_weighted(tour_best_ant, elitist_ants);
           
           if (global_best_tour > distanceForeing){
           	global_best_tour = distanceForeing;
           	/*write_best_global_tour(global_best_tour);*/
           }
           
        } else {
        
           if(distanceOwn < global_best_tour){
              global_best_tour = distanceOwn;
            /*  write_best_global_tour(global_best_tour);*/
           }
        }
     }
  }
  
}


/*******************************************************************/
/************************** Asyncronous ****************************/

/* UNA VEZ PROBADO EL MÉTODO SINCRONO PUEDES IMPLEMENTAR EL ASÍNCRONO, QUE ES MÁS COMPLICADO.
 EN ESTE CASO SE TRATA DE QUE CADA VEZ QUE UNA COLONIA DETECTA UNA SOLUCIÓN MEJOR (ESTO ESTÁ 
 EN LA FUNCIÓN "update_statistics") SE ENVÍA LA RUTA A TODAS LAS DEMÁS COLONIAS, PERO DE FORMA 
 ASÍNCRONA (MPI_Isend). ASI EL PROCESO NO SE QUEDA BLOQUEADO A LA ESPERA DE QUE LOS DEMÁS LO 
 RECIBAN.
 
 TODOS LOS PROCESOS HAN DE INICIALIZAR AL COMIENZO UN BUFFER PREPARADO PARA RECIBIR COMUNICACIONES
 EN CUALQUIER MOMENTO (MPI_Irecv) QUE TAMPOCO BLOQUEA LOS PROCESOS.
 
 ANTES DE TERMINAR CADA ITERACIÓN (EN LA FUNCIÓN termination_condition) SE PUEDE CHEQUEAR SI
 SE HAN RECIBIDO NUEVAS MEJORES RUTAS DEL RESTO DE COLONIAS. SI ES ASÍ SE ACTUALIZA LA MATRIZ DE 
 FEROMONAS, SI NO SE HAN RECIBIDO SEGUIMOS ADELANTE CON LA EJECUCIÓN.
*/

/*
 Routine to start recv. of communications from Colonies
 
 AQUÍ SE PREPARA EL BUFFER DE RECEPCIÓN ASÍNCRONA (MPI_Irecv)
 
*/
/*void startCommColonies ( void )
{*/
    /* Prepare to receive the best solutions found in
     the rest of the colonies */
 /*  MPI_Irecv(XXXXX);
}*/


/*
 Routine to send the best solution found to the rest of Colonies
 
 AQUÍ SE ENVIARÁ CON MPI_Isend LA MEJOR SOLUCIÓN ENCONTRADA AL RESTO DE COLONIAS
*/
/*void sendBestSolutionToColonies ( void )
{
 MPI_Isend
} */

/*
 Routine to check if there are pending messages from Colonies
 
 AQUÍ HABRÁ UN LAZO WHILE QUE ESCUCHARÁ LOS MENSAJES (PORQUE PUEDE HABER VARIOS 
 PENDIENTES) CON UN MPI_Test. CUANDO YA NO QUEDEN MENSAJES PENDIENTES SEGUIRÁ 
 ADELANTE CON LA EJECUCIÓN.
 
*/
/*void listenColonies()
{

MPI_Test
 
}*/

/*************************************************************************/
/********************* Update pheromone **********************************/

/*
 Routines to update pheromone using the foreigner tour
*/
void foreign_solution_update_pheromone( long int *ftour )
/*
 FUNCTION:      reinforces edges used in foreigner solution
 INPUT:         pointer to foreigner tour that updates the pheromone trail
 OUTPUT:        none
 (SIDE)EFFECTS: pheromones of arcs in foreigner tour are increased
 */
{
    long int i, j, h, ftour_length;
    double   d_tau;
    
    TRACE ( printf("global pheromone update with foreigner solutions\n"); );
    
    ftour_length = compute_tour_length( ftour );
    d_tau = 1.0 / (double) ftour_length;
    for ( i = 0 ; i < n ; i++ ) {
        j = ftour[i];
        h = ftour[i+1];
        pheromone[j][h] += d_tau;
        pheromone[h][j] = pheromone[j][h];
    }
}

void foreign_solution_update_pheromone_weighted( long int *ftour, long int weight )
/*
 FUNCTION:      reinforces edges used in foreigner solution with weight
 INPUT:         pointer to foreigner tour that updates the pheromone trail and its weight
 OUTPUT:        none
 (SIDE)EFFECTS: pheromones of arcs in foreigner tour are increased
 */
{
    long int      i, j, h, ftour_length;
    double        d_tau;
    
    TRACE ( printf("global pheromone update with foreigner solutions weighted\n"); );
    
    ftour_length = compute_tour_length( ftour );
    d_tau = (double) weight / (double) ftour_length;
    for ( i = 0 ; i < n ; i++ ) {
        j = ftour[i];
        h = ftour[i+1];
        pheromone[j][h] += d_tau;
        pheromone[h][j] = pheromone[j][h];
    }       
}


