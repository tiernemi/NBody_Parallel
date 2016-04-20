/*
 * =====================================================================================
 *
 *       Filename:  mpi_utils.cpp
 *
 *    Description:  Source file for mpi_utils.
 *
 *        Version:  1.0
 *        Created:  04/19/2016 12:11:37 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

#include "mpi_utils.hpp"
#include <mpi.h>

int MPIUtils::rank ;
int MPIUtils::nproc ;

/* 
 * ===  MEMBER FUNCTION CLASS : MPIUtils  ======================================
 *         Name:  initMPI
 *    Arguments:  int argc - Number of arguments.
 *                char * argv[] array of cmd args.
 *  Description:  Initialises MPI and sets rank and nproc.
 * =====================================================================================
 */

void MPIUtils::initMPI(int argc, char * argv[]) {
	MPI_Init(&argc, &argv) ;
	int rankl,nprocl ;
	MPI_Comm_rank(MPI_COMM_WORLD, &rankl) ;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocl) ;
	MPIUtils::rank = rankl ;
	MPIUtils::nproc = nprocl ;
}		/* -----  end of member function initMPI  ----- */


/* 
 * ===  MEMBER FUNCTION CLASS : MPIUtils  ======================================
 *         Name:  finaliseMPI
 *  Description:  Finalises MPI.
 * =====================================================================================
 */

void MPIUtils::finaliseMPI() {
	MPI_Finalize() ;
}		/* -----  end of member function finaliseMPI  ----- */
