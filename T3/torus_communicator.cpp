/*
 * =====================================================================================
 *
 *       Filename:  torus_communicator.cpp
 *
 *    Description:  Source file for torus communicator.
 *
 *        Version:  1.0
 *        Created:  04/21/2016 10:32:16 AM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

#include <mpi.h>
#include "torus_communicator.hpp"

MPI_Comm TorusCommunicator::communicator ;
int TorusCommunicator::cartCoords[2] ;
int TorusCommunicator::commDims[2] = {0,0} ;
int TorusCommunicator::neighbourRanks[8] ;

/* 
 * ===  MEMBER FUNCTION CLASS : TorusCommunicator  =====================================
 *         Name:  initCommunicator
 *  Description:  Initialise the cartesian communicator, generates dimensions based on
 *                the number of processes and determines the ranks of the neighbours.
 *                Also enforces periodic boundary conditions.
 * =====================================================================================
 */

void TorusCommunicator::initCommunicator() {
	// Get dims and co-ords as well as cretae communicator. //
	MPI_Dims_create(MPIUtils::nproc,2,commDims) ;
	int periodic[2] = {1,1} ;
	MPI_Cart_create(MPI_COMM_WORLD,2,commDims,periodic,true,&communicator) ;
	MPI_Cart_coords(communicator,MPIUtils::rank,2,cartCoords);
	// Find the neighbouring ranks. //
	MPI_Cart_shift(communicator,1,1,&neighbourRanks[LEFT],&neighbourRanks[RIGHT]) ;
	MPI_Cart_shift(communicator,0,1,&neighbourRanks[BOT],&neighbourRanks[TOP]) ;
	int topLeft[2] ;
	topLeft[0] = cartCoords[0]+1 ;
	topLeft[1] = cartCoords[1]-1 ;
	MPI_Cart_rank(communicator,topLeft,&neighbourRanks[TOPLEFT]);
	int topRight[2] ;
	topRight[0] = cartCoords[0]+1 ;
	topRight[1] = cartCoords[1]+1 ;
	MPI_Cart_rank(communicator,topRight,&neighbourRanks[TOPRIGHT]);
	int botLeft[2] ;
	botLeft[0] = cartCoords[0]-1 ;
	botLeft[1] = cartCoords[1]-1 ;
	MPI_Cart_rank(communicator,botLeft,&neighbourRanks[BOTLEFT]);
	int botRight[2] ;
	botRight[0] = cartCoords[0]-1 ;
	botRight[1] = cartCoords[1]+1 ;
	MPI_Cart_rank(communicator,botRight,&neighbourRanks[BOTRIGHT]);
}

/* 
 * ===  MEMBER FUNCTION CLASS : TorusCommunicator  ======================================
 *         Name:  sendRecvSizes.
 *    Arguments:  int source - The source of the size information.
 *                int dest - The destination of the size information.
 *                int & numSent - Number of data points sent to dest.
 *                int * numRecv - Number of data points to recieve from dest.
 *                MPI_Request & reqPS - Request for sending.
 *                MPI_Request & reqPR - Request for recieving.
 *  Description:  Send and recieve the size of the data sends required between source and
 *                dest.
 * =====================================================================================
 */

void TorusCommunicator::sendRecvSizes(int source, int dest, int & numSent, int & numRecv, MPI_Request & reqPS, MPI_Request & reqPR) {
	MPI_Isend(&numSent,1,MPI_INT,neighbourRanks[source],source,communicator,&reqPS);
	MPI_Request_free(&reqPS) ;
	MPI_Irecv(&numRecv,1,MPI_INT,neighbourRanks[dest],source,communicator,&reqPR) ;
}		/* -----  end of member function   ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : TorusCommunicator  ======================================
 *         Name:  getCommunicator
 *      Returns:  Cartesian communicator.
 *  Description:  Returns the communicator.
 * =====================================================================================
 */

MPI_Comm & TorusCommunicator::getCommunicator() {
	return communicator ;
}		/* -----  end of member function getCommunicator  ----- */


/* 
 * ===  MEMBER FUNCTION CLASS : TorusCommunicator  ======================================
 *         Name:  getCoords.
 *      Returns:  The co-ordinates array.
 *  Description:  Returns array of cartesian co-ordinates of system.
 * =====================================================================================
 */

int * TorusCommunicator::getCoords() {
	return cartCoords ;
}		/* -----  end of member function getCommunicator  ----- */


/* 
 * ===  MEMBER FUNCTION CLASS : TorusCommunicator  ======================================
 *         Name:  getNeighbours
 *      Returns:  An array of nearest neighbour ranks.
 * =====================================================================================
 */

int * TorusCommunicator::getNeighbours() {
	return neighbourRanks ;
}		/* -----  end of member function   ----- */


/* 
 * ===  MEMBER FUNCTION CLASS : TorusCommunicator  ======================================
 *         Name:  getCommDims
 *      Returns:  An array of the communicator dimensions.
 * =====================================================================================
 */

int * TorusCommunicator::getCommDims() {
	return commDims ;
}		/* -----  end of member function   ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : TorusCommunicator  ======================================
 *         Name:  destroyCommunicator
 *  Description:  Frees memory associated with communicator.
 * =====================================================================================
 */

void TorusCommunicator::destroyCommunicator() {
	MPI_Comm_free(&communicator) ;
}		/* -----  end of member function function  ----- */
