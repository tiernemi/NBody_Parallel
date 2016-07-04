#ifndef MPI_UTILS_HPP_QV2JD7K1
#define MPI_UTILS_HPP_QV2JD7K1

/*
 * =====================================================================================
 *
 *       Filename:  mpi_utils.hpp
 *
 *    Description:  Static helper class for MPI utilities
 *
 *        Version:  1.0
 *        Created:  04/19/2016 12:16:17 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

/* 
 * ===  CLASS  =========================================================================
 *         Name:  MPIUtils
 *       Fields:  int rank - Rank of process.
 *                int nproc - Number of processes.
 *  Description:  Barebones helper class for MPI functions.
 * =====================================================================================
 */

class MPIUtils {
 public:
	 static void initMPI(int argc, char * argv[]) ;
	 static void finaliseMPI() ;
 public:
	static int rank ;
	static int nproc ;
} ;		/* -----  end of class MPIUtils  ----- */


#endif /* end of include guard: MPI_UTILS_HPP_QV2JD7K1 */
