#ifndef TORUS_COMMUNICATOR_HPP_HG0MKCBV
#define TORUS_COMMUNICATOR_HPP_HG0MKCBV

/*
 * =====================================================================================
 *
 *       Filename:  torus_communicator.hpp
 *
 *    Description:  Communicator policy for N-Body system.
 *
 *        Version:  1.0
 *        Created:  04/19/2016 12:36:45 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

#include "mpi.h"
#include <vector>
#include "mpi_utils.hpp"

/* 
 * ===  CLASS  =========================================================================
 *         Name:  TorusCommunicator
 *       Fields:  
 *  Description:  
 * =====================================================================================
 */

class TorusCommunicator {
 public:
	static void initCommunicator() ; 
	static void destroyCommunicator() ;
	static MPI_Comm & getCommunicator() ;
	static int * getCoords() ;
	static int * getNeighbours() ;
	static int * getCommDims() ;
	template <typename DataType>
	static inline void communicateBoundaries(std::vector<std::vector<DataType>> &, std::vector<std::vector<DataType>> &) ;
	template <typename DataType>
	static inline void tradeParticles(std::vector<std::vector<DataType>> &, std::vector<std::vector<DataType>> &) ;
	// Helper enum. //
	enum NEIGHBOURS {
		LEFT,
		RIGHT,
		TOP,
		BOT,
		TOPLEFT,
		TOPRIGHT,
		BOTLEFT,
		BOTRIGHT
	} ;
	static NEIGHBOURS neighEnum ;
 private:
	void sendSizes(int, int, int *, int *, MPI_Request &, MPI_Request &) ;
	template <typename DataType>
	void sendData(int, int, std::vector<std::vector<DataType>> &, std::vector<std::vector<DataType>> &, MPI_Request &, MPI_Request &) ;
 private:
	static MPI_Comm communicator ;
	static int cartCoords[2] ;
	static int commDims[2] ;
	static int neighbourRanks[8] ;
	/* data */
} ;		/* -----  end of class TorusCommunicator  ----- */

// TEMPLATED MEMBER FUNCTIONS //

/* 
 * ===  MEMBER FUNCTION CLASS : TorusCommunicator  ======================================
 *         Name:  commuincateBoundaries
 *    Arguments:  int direction - Direction in which to send particle data.
 *                Datatype * sendbuf - Buffer used for sending data.
 *  Description:  Sends particles in send buffer in specified direction.
 * =====================================================================================
 */

template <typename DataType>
void inline TorusCommunicator::communicateBoundaries(std::vector<std::vector<DataType>> & sendBuffers, std::vector<std::vector<DataType>> & recBuffers) {
	MPI_Request reqS[8] ;
	MPI_Request reqR[8] ;
	MPI_Request reqPS[8] ;
	MPI_Request reqPR[8] ;

	int numElementsSent[8] ;
	int numElementsRecv[8] ;
	for (int i = 0 ; i < 8 ; ++i) {
		numElementsSent[i] = sendBuffers[i].size() ;
	}
	// Send data sizes left. //
	MPI_Isend(&numElementsSent[LEFT],1,MPI_INT,neighbourRanks[LEFT],LEFT,communicator,&reqPS[LEFT]);
	MPI_Request_free(&reqPS[LEFT]) ;
	MPI_Irecv(&numElementsRecv[RIGHT],1,MPI_INT,neighbourRanks[RIGHT],LEFT,communicator,&reqPR[RIGHT]) ;
	// Send data sizes left. //
	MPI_Isend(&numElementsSent[RIGHT],1,MPI_INT,neighbourRanks[RIGHT],RIGHT,communicator,&reqPS[RIGHT]);
	MPI_Request_free(&reqPS[RIGHT]) ;
	MPI_Irecv(&numElementsRecv[LEFT],1,MPI_INT,neighbourRanks[LEFT],RIGHT,communicator,&reqPR[LEFT]) ;
	// Send data sizes left. //
	MPI_Isend(&numElementsSent[TOP],1,MPI_INT,neighbourRanks[TOP],TOP,communicator,&reqPS[TOP]);
	MPI_Request_free(&reqPS[TOP]) ;
	MPI_Irecv(&numElementsRecv[BOT],1,MPI_INT,neighbourRanks[BOT],TOP,communicator,&reqPR[BOT]) ;
	// Send data sizes left. //
	MPI_Isend(&numElementsSent[BOT],1,MPI_INT,neighbourRanks[BOT],BOT,communicator,&reqPS[BOT]);
	MPI_Request_free(&reqPS[BOT]) ;
	MPI_Irecv(&numElementsRecv[TOP],1,MPI_INT,neighbourRanks[TOP],BOT,communicator,&reqPR[TOP]) ;

	MPI_Waitall(8,reqPR,MPI_STATUS_IGNORE) ;

	// Send Left boundaries. //
	if (numElementsSent[LEFT] != 0) {
		MPI_Isend(&sendBuffers[LEFT][0],numElementsSent[LEFT],MPI_DOUBLE,neighbourRanks[LEFT],LEFT,communicator,&reqS[LEFT]);
		MPI_Request_free(&reqS[LEFT]) ;
	}
	if (numElementsRecv[RIGHT] != 0) {
		recBuffers[RIGHT].resize(numElementsRecv[RIGHT]) ;
		MPI_Irecv(&recBuffers[RIGHT][0],numElementsRecv[RIGHT],MPI_DOUBLE,neighbourRanks[RIGHT],LEFT,communicator,&reqR[RIGHT]) ;
	}
	// Send RIGHT boundaries. //
	if (numElementsSent[RIGHT] != 0) {
		MPI_Isend(&sendBuffers[RIGHT][0],numElementsSent[RIGHT],MPI_DOUBLE,neighbourRanks[RIGHT],RIGHT,communicator,&reqS[RIGHT]);
		MPI_Request_free(&reqS[RIGHT]) ;
	}
	if (numElementsRecv[LEFT] != 0) {
		recBuffers[LEFT].resize(numElementsRecv[LEFT]) ;
		MPI_Irecv(&recBuffers[LEFT][0],numElementsRecv[LEFT],MPI_DOUBLE,neighbourRanks[LEFT],RIGHT,communicator,&reqR[LEFT]) ;
	}
	// Send BOT boundaries. //
	if (numElementsSent[BOT] != 0) {
		MPI_Isend(&sendBuffers[BOT][0],numElementsSent[BOT],MPI_DOUBLE,neighbourRanks[BOT],BOT,communicator,&reqS[BOT]);
		MPI_Request_free(&reqS[BOT]) ;
	}
	if (numElementsRecv[TOP] != 0) {
		recBuffers[TOP].resize(numElementsRecv[TOP]) ;
		MPI_Irecv(&recBuffers[TOP][0],numElementsRecv[TOP],MPI_DOUBLE,neighbourRanks[TOP],BOT,communicator,&reqR[TOP]) ;
	}
	// Send TOP boundaries. //
	if (numElementsSent[TOP] != 0) {
		MPI_Isend(&sendBuffers[TOP][0],numElementsSent[TOP],MPI_DOUBLE,neighbourRanks[TOP],TOP,communicator,&reqS[TOP]);
		MPI_Request_free(&reqS[TOP]) ;
	}
	if (numElementsRecv[BOT] != 0) {
		recBuffers[BOT].resize(numElementsRecv[BOT]) ;
		MPI_Irecv(&recBuffers[BOT][0],numElementsRecv[BOT],MPI_DOUBLE,neighbourRanks[BOT],TOP,communicator,&reqR[BOT]) ;
	}
	// Send TOPLEFT boundaries. //
	if (numElementsSent[TOPLEFT] != 0) {
		MPI_Isend(&sendBuffers[TOPLEFT][0],numElementsSent[TOPLEFT],MPI_DOUBLE,neighbourRanks[TOPLEFT],TOPLEFT,communicator,&reqS[TOPLEFT]);
		MPI_Request_free(&reqS[TOPLEFT]) ;
	}
	if (numElementsRecv[BOTRIGHT] != 0) {
		recBuffers[BOTRIGHT].resize(numElementsRecv[BOTRIGHT]) ;
		MPI_Irecv(&recBuffers[BOTRIGHT][0],numElementsRecv[BOTRIGHT],MPI_DOUBLE,neighbourRanks[BOTRIGHT],TOPLEFT,communicator,&reqR[BOTRIGHT]) ;
	}
	// Send TOPRIGHT boundaries. //
	if (numElementsSent[TOPRIGHT] != 0) {
		MPI_Isend(&sendBuffers[TOPRIGHT][0],numElementsSent[TOPRIGHT],MPI_DOUBLE,neighbourRanks[TOPRIGHT],TOPRIGHT,communicator,&reqS[TOPRIGHT]);
		MPI_Request_free(&reqS[TOPRIGHT]) ;
	}
	if (numElementsRecv[BOTLEFT] != 0) {
		recBuffers[BOTLEFT].resize(numElementsRecv[BOTLEFT]) ;
		MPI_Irecv(&recBuffers[BOTLEFT][0],numElementsRecv[BOTLEFT],MPI_DOUBLE,neighbourRanks[BOTLEFT],TOPRIGHT,communicator,&reqR[BOTLEFT]) ;
	}
	// Send BOTLEFT boundaries. //
	if (numElementsSent[BOTLEFT] != 0) {
		MPI_Isend(&sendBuffers[BOTLEFT][0],numElementsSent[BOTLEFT],MPI_DOUBLE,neighbourRanks[BOTLEFT],BOTLEFT,communicator,&reqS[BOTLEFT]);
		MPI_Request_free(&reqS[BOTLEFT]) ;
	}
	if (numElementsRecv[TOPRIGHT] != 0) {
		recBuffers[TOPRIGHT].resize(numElementsRecv[TOPRIGHT]) ;
		MPI_Irecv(&recBuffers[TOPRIGHT][0],numElementsRecv[TOPRIGHT],MPI_DOUBLE,neighbourRanks[TOPRIGHT],BOTLEFT,communicator,&reqR[TOPRIGHT]) ;
	}
	// Send BOTRIGHT boundaries. //
	if (numElementsSent[BOTRIGHT] != 0) {
		MPI_Isend(&sendBuffers[BOTRIGHT][0],numElementsSent[BOTRIGHT],MPI_DOUBLE,neighbourRanks[BOTRIGHT],BOTRIGHT,communicator,&reqS[BOTRIGHT]);
		MPI_Request_free(&reqS[BOTRIGHT]) ;
	}
	if (numElementsRecv[TOPLEFT] != 0) {
		recBuffers[TOPLEFT].resize(numElementsRecv[TOPLEFT]) ;
		MPI_Irecv(&recBuffers[TOPLEFT][0],numElementsRecv[TOPLEFT],MPI_DOUBLE,neighbourRanks[TOPLEFT],BOTRIGHT,communicator,&reqR[TOPLEFT]) ;
	}	
	MPI_Waitall(8,reqR,MPI_STATUS_IGNORE) ;
}

/*  
 
 *===  MEMBER FUNCTION CLASS : TorusCommunicator  ======================================
 *         Name:  recvBoundaryParticle
 *    Arguments:  int direction - Direction from which to recieve particle data.
 *                Datatype * sendbuf - Buffer used for receiving data.
 *  Description:  Recieves particles into buffer from specified direction.
 * =====================================================================================
 

template <typename DataType>
void inline TorusCommunicator::recvBoundaryParticles(int direction, DataType * sendBuffer) {
}
*/


#endif /* end of include guard: TORUS_COMMUNICATOR_HPP_HG0MKCBV */

