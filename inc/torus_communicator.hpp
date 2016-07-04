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
 *       Fields:  MPI_Comm communicator - The cartesian communicator for system.
 *                int cartCoords[2] - The cartesian co-ordinates of system.
 *                int commDims[2] - The communicator dimensions.
 *                int neighbourRanks[8] - The ranks of the neighbouring processes.
 *  Description:  Communicator used to send/recieve buffer data to and from 8 nearest
 *                neighbours.
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
	static void sendRecvSizes(int, int, int &, int &, MPI_Request &, MPI_Request &) ;
	template <typename DataType>
	static void sendRecvData(int, int, std::vector<DataType> &, std::vector<DataType> &, int &, int &, MPI_Request &, MPI_Request &) ;
 private:
	static MPI_Comm communicator ;
	static int cartCoords[2] ;
	static int commDims[2] ;
	static int neighbourRanks[8] ;
	/* data */
} ;		/* -----  end of class TorusCommunicator  ----- */

// TEMPLATED MEMBER FUNCTIONS //

/* 
 * ===  MEMBER FUNCTION CLASS : TorusCommunicator  ========================================
 *         Name:  commuincateBoundaries
 *    Arguments:  int direction - Direction in which to send particle data.
 *                Datatype * sendbuf - Buffer used for sending data.
 *  Description:  Send corresponding send buffers to 8 nearest neighbours and then
 *                recieves data from these 8 nearest neighbours.
 * ============================================================================================
 */

template <typename DataType>
void inline TorusCommunicator::communicateBoundaries(std::vector<std::vector<DataType>> & sendBuffers, std::vector<std::vector<DataType>> & recvBuffers) {
	MPI_Request reqS[8] ;
	MPI_Request reqR[8] ;
	MPI_Request reqPS[8] ;
	MPI_Request reqPR[8] ;

	// Initialise request objects.
	for (int i = 0 ; i < 8 ; ++i) {	
		reqS[i] = MPI_REQUEST_NULL ;
		reqR[i] = MPI_REQUEST_NULL ;
		reqPS[i] = MPI_REQUEST_NULL ;
		reqPR[i] = MPI_REQUEST_NULL ;
	}

	int numElementsSent[8] ;
	int numElementsRecv[8] ;
	for (int i = 0 ; i < 8 ; ++i) {
		numElementsSent[i] = sendBuffers[i].size() ;
	}

	// Send and recieve the size of the data buffers for each neighbour. //
	sendRecvSizes(LEFT,RIGHT,numElementsSent[LEFT],numElementsRecv[RIGHT],reqPS[LEFT],reqPR[RIGHT]) ;
	sendRecvSizes(RIGHT,LEFT,numElementsSent[RIGHT],numElementsRecv[LEFT],reqPS[RIGHT],reqPR[LEFT]) ;
	sendRecvSizes(BOT,TOP,numElementsSent[BOT],numElementsRecv[TOP],reqPS[BOT],reqPR[TOP]) ;
	sendRecvSizes(TOP,BOT,numElementsSent[TOP],numElementsRecv[BOT],reqPS[TOP],reqPR[BOT]) ;
	sendRecvSizes(TOPLEFT,BOTRIGHT,numElementsSent[TOPLEFT],numElementsRecv[BOTRIGHT],reqPS[TOPLEFT],reqPR[BOTRIGHT]) ;
	sendRecvSizes(TOPRIGHT,BOTLEFT,numElementsSent[TOPRIGHT],numElementsRecv[BOTLEFT],reqPS[TOPRIGHT],reqPR[BOTLEFT]) ;
	sendRecvSizes(BOTLEFT,TOPRIGHT,numElementsSent[BOTLEFT],numElementsRecv[TOPRIGHT],reqPS[BOTLEFT],reqPR[TOPRIGHT]) ;
	sendRecvSizes(BOTRIGHT,TOPLEFT,numElementsSent[BOTRIGHT],numElementsRecv[TOPLEFT],reqPS[BOTRIGHT],reqPR[TOPLEFT]) ;
	MPI_Waitall(8,reqPR,MPI_STATUS_IGNORE) ;

	// Send and recieve the data buffers for each neighbour. //
	sendRecvData(LEFT,RIGHT,sendBuffers[LEFT],recvBuffers[RIGHT],numElementsSent[LEFT],numElementsRecv[RIGHT],reqS[LEFT],reqR[RIGHT]) ;
	sendRecvData(RIGHT,LEFT,sendBuffers[RIGHT],recvBuffers[LEFT],numElementsSent[RIGHT],numElementsRecv[LEFT],reqS[RIGHT],reqR[LEFT]) ;
	sendRecvData(BOT,TOP,sendBuffers[BOT],recvBuffers[TOP],numElementsSent[BOT],numElementsRecv[TOP],reqS[BOT],reqR[TOP]) ;
	sendRecvData(TOP,BOT,sendBuffers[TOP],recvBuffers[BOT],numElementsSent[TOP],numElementsRecv[BOT],reqS[TOP],reqR[BOT]) ;
	sendRecvData(TOPLEFT,BOTRIGHT,sendBuffers[TOPLEFT],recvBuffers[BOTRIGHT],numElementsSent[TOPLEFT],numElementsRecv[BOTRIGHT],reqS[TOPLEFT],reqR[BOTRIGHT]) ;
	sendRecvData(TOPRIGHT,BOTLEFT,sendBuffers[TOPRIGHT],recvBuffers[BOTLEFT],numElementsSent[TOPRIGHT],numElementsRecv[BOTLEFT],reqS[TOPRIGHT],reqR[BOTLEFT]) ;
	sendRecvData(BOTLEFT,TOPRIGHT,sendBuffers[BOTLEFT],recvBuffers[TOPRIGHT],numElementsSent[BOTLEFT],numElementsRecv[TOPRIGHT],reqS[BOTLEFT],reqR[TOPRIGHT]) ;
	sendRecvData(BOTRIGHT,TOPLEFT,sendBuffers[BOTRIGHT],recvBuffers[TOPLEFT],numElementsSent[BOTRIGHT],numElementsRecv[TOPLEFT],reqS[BOTRIGHT],reqR[TOPLEFT]) ;

	MPI_Waitall(8,reqR,MPI_STATUS_IGNORE) ;
	MPI_Barrier(MPI_COMM_WORLD) ;
}

 /* 
 * ===  MEMBER FUNCTION CLASS : TorusCommunicator  ==============================================
 *         Name:  sendRecvData.
 *    Arguments:  int source - The source of the size information.
 *                int dest - The destination of the size information.
 *                std::vector<Datatype> & sendData - Transition and/or ghost data to be sent.
 *                std::vector<Datatype> & recvData - Transition and/or ghost data to be recieved.
 *                int & numSent - Number of data points sent to dest.
 *                int * numRecv - Number of data points to recieve from dest.
 *                MPI_Request & reqPS - Request for sending.
 *                MPI_Request & reqPR - Request for recieving.
 *  Description:  Send and recieve the transition and ghost data between process source and dest.
 * ==============================================================================================
 */


template <typename DataType>
void TorusCommunicator::sendRecvData(int source, int dest, std::vector<DataType> & sendData, 
		           std::vector<DataType> & recvData, int & numSent, int & numRecv, MPI_Request & reqS, MPI_Request & reqR) {
	// Send data to dest. //
	MPI_Isend(&sendData[0],numSent,MPI_DOUBLE,neighbourRanks[source],source,communicator,&reqS);
	MPI_Request_free(&reqS) ;
	recvData.resize(numRecv) ;
	// Recieve data from dest. //
	MPI_Irecv(&recvData[0],numRecv,MPI_DOUBLE,neighbourRanks[dest],source,communicator,&reqR) ;
}		/* -----  end of member function function  ----- */


#endif /* end of include guard: TORUS_COMMUNICATOR_HPP_HG0MKCBV */

