#ifndef SYSTEM_HPP_PSH64OVQ
#define SYSTEM_HPP_PSH64OVQ

/*
 * =====================================================================================
 *
 *       Filename:  system.hpp
 *
 *    Description:  Header file for system class. The system contains N particles and
 *                  evaluates them based on an evaluation policy and updates based
 *                  on an update policy.
 *
 *        Version:  1.0
 *        Created:  04/12/2016 01:20:59 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include "mpi_utils.hpp"
#include <algorithm>

/* 
 * ===  CLASS  =========================================================================
 *         Name:  System
 *       Fields:  std::vector<DataType> positions - Array storing position data.
 *	              std::vector<DataType> velocities - Array storing velocity data.
 *	              std::vector<DataType> forces - Array storing force data.
 *  Description:  System of particles subject to a particular potential specified by the
 *                policy Potential and a particular integration scheme specified by the
 *                policy Integrator. All data associated with the particles is stored
 *                in contiguous arrays within in the system class. 
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
class System {
 public:
	System(const DataType & cellLength) ;
	virtual ~System() ;
	void randomInitialise(int, int) ;
	void addParticle(const std::pair<DataType,DataType> &, const std::pair<DataType,DataType> &) ;
	void simulate(int) ;
	void simulate(int, std::ostream &) ;
	void simulateEnergies(int, std::ostream &, std::ostream &) ;
	void simulateEnergies(int, std::ostream &) ;
	int numParticles() const ;
	void printSystem(std::ostream &, int iter) const ;
	void printEnergies(std::ostream &, int iter) const ;
 private:
	void updatePositions() ;
	void updateVelocities() ;
	void updateVelocitiesEnergies() ;
	void confineParticles() ;
	void updateForces() ;
	void updateForcesEnergies() ;
	void appendInteractionsToBuffers() ;
	void appendTransitionsToBuffers() ;
	void communicateInteractionsTransitions() ;
	void clearBuffers() ;
	void unpackData() ;
 private:
	const DataType totCellLength ;
	const DataType totCellLengthInv ;
	DataType cellLengthX ;
	DataType cellLengthY ;
	DataType lDxZone ;
	DataType rDxZone ;
	DataType tDyZone ;
	DataType bDyZone ;
	DataType lGlobalXOffset ;
	DataType bGlobalYOffset ;
	DataType rGlobalXOffset ;
	DataType tGlobalYOffset ;
	std::vector<DataType> positions ;
	std::vector<DataType> velocities ;
	std::vector<DataType> forces ;
	std::vector<DataType> boundaryGhostPositions ;
	std::vector<std::vector<DataType>> sendBuffers ;
	std::vector<std::vector<DataType>> recvBuffers ;
	DataType totEnergy ;
	DataType initialEnergy ;
} ;		/* -----  end of class System  ----- */

// TEMPLATED MEMBER FUNCTIONS //

/* 
 * ===  MEMBER FUNCTION CLASS : System  ======================================
 *         Name:  System
 *    Arguments:  
 *      Returns:  
 *  Description:  
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
System<DataType,Integrator,Potential,Communicator>::System(const DataType & totCellLength) : totCellLength{totCellLength}, totCellLengthInv{1./totCellLength} {
	Communicator::initCommunicator() ;
	cellLengthY = totCellLength/double(Communicator::getCommDims()[0]) ;
	cellLengthX = totCellLength/double(Communicator::getCommDims()[1]) ;
	bGlobalYOffset = cellLengthY*Communicator::getCoords()[0] ;
	lGlobalXOffset = cellLengthX*Communicator::getCoords()[1] ;
	rGlobalXOffset = lGlobalXOffset + cellLengthX ;
	tGlobalYOffset = bGlobalYOffset + cellLengthY ;
	totEnergy = 0 ;
	sendBuffers.resize(8) ;
	recvBuffers.resize(8) ;
	lDxZone = lGlobalXOffset + Potential::cutOffDist ;
	rDxZone = lGlobalXOffset + cellLengthX - Potential::cutOffDist ;
	bDyZone = bGlobalYOffset + Potential::cutOffDist ;
	tDyZone = bGlobalYOffset + cellLengthY - Potential::cutOffDist ;
}		/* -----  end of member function System  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ===============================================
 *         Name: randomInitialise(int seed)
 *    Arguments: int seed - Seed for rng.
 *               int numParts - Number of particles.
 *  Description: Randomly distributes particles in box given a seed. All particles must
 *               be further than 0.8 apart.
 /===================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
void System<DataType,Integrator,Potential,Communicator>::randomInitialise(int seed, int numParts) {
	std::mt19937 ranGen ;
	ranGen.seed(seed) ;
	std::uniform_real_distribution<DataType> posDistX(lGlobalXOffset+0.3,rGlobalXOffset-0.3) ;
	std::uniform_real_distribution<DataType> posDistY(bGlobalYOffset+0.3,tGlobalYOffset-0.3) ;
	std::normal_distribution<DataType> velDist(0,1) ;
	for (int i = 0 ; i < numParts ; ++i) {
		DataType newPartX = posDistX(ranGen) ;
		DataType newPartY = posDistY(ranGen) ;
		DataType newVelX = velDist(ranGen) ;
		DataType newVelY = velDist(ranGen) ;
		bool tooClose = false ;
		// Check if the particle is too close to another particle. //
		for (unsigned int j = 0 ; j < positions.size() ; j+=2) {
			DataType x_comp = positions[j]-newPartX ;
			x_comp -= totCellLength*std::round(x_comp*totCellLengthInv) ;
			DataType y_comp = positions[j+1]-newPartY ;
			y_comp -= totCellLength*std::round(y_comp*totCellLengthInv) ;
			if (std::sqrt(x_comp*x_comp + y_comp*y_comp) < 0.8) {
				tooClose = true ;
				break ;
			}
		}
		if (tooClose) {
			// Try again.
			--i ;
		} else {
			// Add a new particle. //
			addParticle(std::pair<DataType,DataType>(newPartX,newPartY),
					std::pair<DataType,DataType>(newVelX,newVelY)) ;
		}
	}
}

/* 
 * ===  MEMBER FUNCTION CLASS : System  =================================================
 *         Name:  System
 *    Arguments:  const std::pair<Datatype,Datatype> & position - Position of new particle.
 *                const std::pair<Datatype,Datatype> & velocity - Velocity of new particle.
 *  Description:  Add a new particle to the system.
 * ======================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
void System<DataType,Integrator,Potential,Communicator>::addParticle(const std::pair<DataType,DataType> & position, const std::pair<DataType,DataType> & velocity) {
	positions.push_back(position.first) ;
	positions.push_back(position.second) ;
	velocities.push_back(velocity.first) ;
	velocities.push_back(velocity.second) ;
	forces.push_back(0) ;
	forces.push_back(0) ;
}		/* -----  end of member function System  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  simulate
 *    Arguments:  int numIters - Number of iterations to run simulation.
 *  Description:  Simulate the system for a set number of iterations.
 * ====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
void System<DataType,Integrator,Potential,Communicator>::simulate(int numIters) {
	// Initialise the forces. //
	updateForces() ;
	// Prime the velocity based on the integration policy. //
	for (unsigned int i = 0 ; i < positions.size() ; i+=2) {
		Integrator::initialise(&positions[i], &velocities[i], &forces[i]) ;
	}
	// Update velocities and positions based on the integration policy. //
	for (int i = 0; i < numIters ; ++i) {
		updateVelocities() ;
		updatePositions() ;
		updateForces() ;
		if (i % 1000 == 0) {
			confineParticles() ;
		}
	}
	confineParticles() ;
}		/* -----  end of member function simulate  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  simulate
 *    Arguments:  int numIters - Number of iterations to run simulation.
 *                std::ostream & ouput - Output stream to print positions to.
 *  Description:  Simulate the system for a set number of iterations. Prints positions.
 * ====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
void System<DataType,Integrator,Potential,Communicator>::simulate(int numIters, std::ostream & output) {
	// Prime the velocity based on the integration policy. //
	for (unsigned int i = 0 ; i < positions.size() ; i+=2) {
		Integrator::initialise(&positions[i], &velocities[i], &forces[i]) ;
	}
	// Initialise the forces. //
	clearBuffers() ;
	appendTransitionsToBuffers() ;
	std::cout << 0 << std::endl;
	appendInteractionsToBuffers() ;
	communicateInteractionsTransitions() ;
	updateForces() ;
	// Update velocities and positions based on the integration policy. //
	for (int i = 0; i < numIters ; ++i) {
		updateVelocities() ;
		updatePositions() ;
		clearBuffers() ;
		appendTransitionsToBuffers() ;
		appendInteractionsToBuffers() ;
		communicateInteractionsTransitions() ;
		updateForces() ;
		printSystem(output,i) ;
		output << std::endl << std::endl ;
	}
}		/* -----  end of member function simulate  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  simulateEnergies
 *    Arguments:  int numIters - Number of iterations to run simulation.
 *                std::ostream & ouputPos - Output stream to print positions to.
 *                std::ostream & ouputEn - Output stream to print energies to.
 *  Description:  Simulate the system for a set number of iterations. Prints positions
 *                and energies.
 * ====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
void System<DataType,Integrator,Potential,Communicator>::simulateEnergies(int numIters, std::ostream & outputPos, std::ostream & outputEn) {
	// Initialise the forces. //
	initialEnergy = 0 ;
	totEnergy = 0 ;
	updateForcesEnergies() ;
	// Prime the velocity based on the integration policy. //
	for (unsigned int i = 0 ; i < positions.size() ; i+=2) {
		Integrator::initialise(&positions[i], &velocities[i], &forces[i]) ;
	}

	// Get initial energy of system and advance 1 time step. //
	updateVelocitiesEnergies() ;
	initialEnergy = totEnergy ;
	updatePositions() ;
	printEnergies(outputEn,0) ;
	totEnergy = 0 ;
	updateForcesEnergies() ;
	outputPos << std::endl << std::endl ;
	confineParticles() ;

	// Update velocities and positions based on the integration policy. //
	for (int i = 1; i < numIters ; ++i) {
		updateVelocitiesEnergies() ;
		updatePositions() ;
		printEnergies(outputEn,i) ;
		totEnergy = 0 ;
		updateForcesEnergies() ;
		printSystem(outputPos,i) ;
		outputPos << std::endl << std::endl ;
		if (i % 1000 == 0) {
			confineParticles() ;
		}
	}
	confineParticles() ;
}		/* -----  end of member function simulate  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  simulateEnergies
 *    Arguments:  int numIters - Number of iterations to run simulation.
 *                std::ostream & ouputEn - Output stream to print energies to.
 *  Description:  Simulate the system for a set number of iterations. Prints energies.
 * ====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
void System<DataType,Integrator,Potential,Communicator>::simulateEnergies(int numIters, std::ostream & outputEn) {
	// Initialise the forces. //
	initialEnergy = 0 ;
	totEnergy = 0 ;
	updateForcesEnergies() ;
	// Prime the velocity based on the integration policy. //
	for (unsigned int i = 0 ; i < positions.size() ; i+=2) {
		Integrator::initialise(&positions[i], &velocities[i], &forces[i]) ;
	}

	// Get initial energy of system and advance 1 time step. //
	updateVelocitiesEnergies() ;
	initialEnergy = totEnergy ;
	updatePositions() ;
	totEnergy = 0 ;
	updateForcesEnergies() ;
	confineParticles() ;

	// Update velocities and positions based on the integration policy. //
	for (int i = 1; i < numIters ; ++i) {
		updateVelocitiesEnergies() ;
		updatePositions() ;
		printEnergies(outputEn,i) ;
		totEnergy = 0 ;
		updateForcesEnergies() ;
		if (i % 1000 == 0) {
			confineParticles() ;
		}
	}
	confineParticles() ;
}		/* -----  end of member function simulate  ----- */


/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  numParticles
 *      Returns:  The number of particles in the system.
 *  Description:  Getter function for the number of particles in the system.
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
int System<DataType,Integrator,Potential,Communicator>::numParticles() const {
	return (positions.size()/2) ;
}		/* -----  end of member function num  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  updatePositions
 *  Description:  Update the positions of the particles in the system.
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
void System<DataType,Integrator,Potential,Communicator>::updatePositions() {
	for (unsigned int i = 0 ; i < positions.size() ; i+=2) {
		Integrator::updatePosition(&positions[i],&velocities[i]) ;
	}
}		/* -----  end of member function updatePositions  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  updateVelocities
 *  Description:  Update the velocities of the particles in the system.
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
void System<DataType,Integrator,Potential,Communicator>::updateVelocities() {	
	for (unsigned int i = 0 ; i < positions.size() ; i+=2) {
		Integrator::updateVelocity(&velocities[i],&forces[i]) ;
	}
}		/* -----  end of member function updateVelocities  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  updateVelocitiesEnergies
 *  Description:  Update the velocities of the particles in the system.
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
void System<DataType,Integrator,Potential,Communicator>::updateVelocitiesEnergies() {	
	for (unsigned int i = 0 ; i < positions.size() ; i+=2) {
		Integrator::updateVelocity(&velocities[i],&forces[i]) ;
		totEnergy += (velocities[i]*velocities[i] + velocities[i+1]*velocities[i+1])/2.0 ;
	}
}		/* -----  end of member function updateVelocities  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  ConfineParticles
 *  Description:  Confines the particles to the box in order to enforce periodic boundaries.
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
void System<DataType,Integrator,Potential,Communicator>::confineParticles() {	
	for (unsigned int i = 0 ; i < positions.size() ; i+=2) {
		positions[i] = std::fmod(positions[i]+totCellLength,totCellLength) ; 
		positions[i+1] = std::fmod(positions[i+1]+totCellLength,totCellLength) ; 
	}
}		/* -----  end of member function function  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ======================================
 *         Name:  updateForces
 *  Description:  Update the forces of the particles in the system.
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
void System<DataType,Integrator,Potential,Communicator>::updateForces() {	
	unsigned int numParts = numParticles() ;
	for (auto i = forces.begin(); i != forces.end(); ++i) {
		*i = 0 ;
	}
	// For each ij interaction calculate force. //
	for (unsigned int i = 0 ; i < numParts*2 ; i+=2) {
		for (unsigned int j = i+2 ; j < numParts*2 ; j+=2) {
			// Get shortest connecting vector in periodic system. //
			DataType x_comp = (positions[j] - positions[i]) ;
			x_comp -= totCellLength*std::round(x_comp*totCellLengthInv) ;
			DataType y_comp = positions[j+1] - positions[i+1] ;
			y_comp -= totCellLength*std::round(y_comp*totCellLengthInv) ;
			// Calculate the force for this vector using the distance. //
			DataType distSqr = x_comp*x_comp + y_comp*y_comp ;
			// Cutoff distance is when r = 2.5 sigma  . //
			if (distSqr < Potential::cutOffDistSqr) {  
				DataType forceMag = Potential::calcForce(distSqr) ;
				DataType xForce = forceMag*x_comp ;
				DataType yForce = forceMag*y_comp ;
				forces[i] += xForce ;
				forces[i+1] += yForce ;
				forces[j] -= xForce ;
				forces[j+1] -= yForce ;
			}
		}
	}
	// Consider the ghost particles. //
	for (unsigned int i = 0 ; i < numParts*2 ; i+=2) {
		for (unsigned int j = 0 ; j < boundaryGhostPositions.size() ; j+=2) {
			// Get shortest connecting vector in periodic system. //
			DataType x_comp = (boundaryGhostPositions[j] - positions[i]) ;
			x_comp -= totCellLength*std::round(x_comp*totCellLengthInv) ;
			DataType y_comp = boundaryGhostPositions[j+1] - positions[i+1] ;
			y_comp -= totCellLength*std::round(y_comp*totCellLengthInv) ;
			// Calculate the force for this vector using the distance. //
			DataType distSqr = x_comp*x_comp + y_comp*y_comp ;
			// Cutoff distance is when r = 2.5 sigma  . //
			if (distSqr < Potential::cutOffDistSqr) {  
				DataType forceMag = Potential::calcForce(distSqr) ;
				DataType xForce = forceMag*x_comp ;
				DataType yForce = forceMag*y_comp ;
				forces[i] += xForce ;
				forces[i+1] += yForce ;
			}
		}
	}
}		/* -----  end of member function updateForces  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ======================================
 *         Name:  updateForcesEnergies
 *  Description:  Update the forces and energies of the particles in the system.
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
void System<DataType,Integrator,Potential,Communicator>::updateForcesEnergies() {	
	unsigned int numParts = numParticles() ;
	for (auto i = forces.begin(); i != forces.end(); ++i) {
		*i = 0 ;
	}
	// For each ij interaction calculate force. //
	for (unsigned int i = 0 ; i < numParts ; ++i) {
		for (unsigned int j = i+1 ; j < numParts ; ++j) {	
			// Get shortest connecting vector in periodic system. //
			DataType x_comp = (positions[2*j] - positions[2*i]) ;
			x_comp -= totCellLength*std::round(x_comp*totCellLengthInv) ;
			DataType y_comp = positions[2*j+1] - positions[2*i+1] ;
			y_comp -= totCellLength*std::round(y_comp*totCellLengthInv) ;
			// Calculate the force for this vector using the distance. //
			DataType distSqr = x_comp*x_comp + y_comp*y_comp ;
			// Cutoff distance is when r = 2.5 sigma  . //
			if (distSqr < Potential::cutOffDistSqr) {
				DataType forceMag = Potential::calcForce(distSqr) ;
				DataType xForce = forceMag*x_comp ;
				DataType yForce = forceMag*y_comp ;
				forces[2*i] += xForce ;
				forces[2*i+1] += yForce ;
				forces[2*j] -= xForce ;
				forces[2*j+1] -= yForce ;
				// Calc energies. //
				totEnergy += Potential::calcEnergy(distSqr) ;
			}
		}
	}
}		/* -----  end of member function updateForces  ----- */


/* 
 * ===  MEMBER FUNCTION CLASS : system  ======================================
 *         Name:  function
 *    Arguments:  
 *      Returns:  
 *  Description:  
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
void System<DataType,Integrator,Potential,Communicator>::clearBuffers() {
	for (int i = 0 ; i < 8 ; ++i) {
		sendBuffers[i].clear() ;
		recvBuffers[i].clear() ;
	}
	boundaryGhostPositions.clear() ;
}		/* -----  end of member function function  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  fillDxZoneBuffers
 *  Description:  Detects if particles are near boundaries and fills the buffers,
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
void System<DataType,Integrator,Potential,Communicator>::appendInteractionsToBuffers() {
	for (unsigned int i = 0 ; i < positions.size() ; i+=2) {
		DataType x_comp = positions[i] ;
		DataType y_comp = positions[i+1] ;
		if (x_comp < lDxZone) {
			sendBuffers[Communicator::LEFT].push_back(positions[i]) ;
			sendBuffers[Communicator::LEFT].push_back(positions[i+1]) ;
			if (y_comp < bDyZone) {
				sendBuffers[Communicator::BOTLEFT].push_back(positions[i]) ;
				sendBuffers[Communicator::BOTLEFT].push_back(positions[i+1]) ;
			} else if (y_comp > tDyZone) {
				sendBuffers[Communicator::TOPLEFT].push_back(positions[i]) ;
				sendBuffers[Communicator::TOPLEFT].push_back(positions[i+1]) ;
			}
		} else if (x_comp > rDxZone) {
			sendBuffers[Communicator::RIGHT].push_back(positions[i]) ;
			sendBuffers[Communicator::RIGHT].push_back(positions[i+1]) ;
			if (y_comp < bDyZone) {
				sendBuffers[Communicator::BOTRIGHT].push_back(positions[i]) ;
				sendBuffers[Communicator::BOTRIGHT].push_back(positions[i+1]) ;
			} else if (y_comp > tDyZone) {
				sendBuffers[Communicator::TOPRIGHT].push_back(positions[i]) ;
				sendBuffers[Communicator::TOPRIGHT].push_back(positions[i+1]) ;
			}
		}
		if (y_comp < bDyZone) {
			sendBuffers[Communicator::BOT].push_back(positions[i]) ;
			sendBuffers[Communicator::BOT].push_back(positions[i+1]) ;
		} else if (y_comp > tDyZone) {
			sendBuffers[Communicator::TOP].push_back(positions[i]) ;
			sendBuffers[Communicator::TOP].push_back(positions[i+1]) ;
		}
	}
}		/* -----  end of member function function  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  fillDxZoneBuffers
 *  Description:  Detects if particles are near boundaries and fills the buffers,
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
void System<DataType,Integrator,Potential,Communicator>::appendTransitionsToBuffers() {
	// Find transitions and add to relevent buffer. //
	for (unsigned int i = 0 ; i < positions.size() ; i+=2) {
		DataType x_comp = positions[i] ;
		DataType y_comp = positions[i+1] ;
		bool hasLeftProcess = false ;
		if (x_comp < lGlobalXOffset) {
			if (x_comp < 0) {
				positions[i] = totCellLength + positions[i] ;
			}
			sendBuffers[Communicator::LEFT].push_back(positions[i]) ;
			sendBuffers[Communicator::LEFT].push_back(positions[i+1]) ;
			sendBuffers[Communicator::LEFT].push_back(velocities[i]) ;
			sendBuffers[Communicator::LEFT].push_back(velocities[i+1]) ;
			hasLeftProcess = true ;
		} else if (x_comp > rGlobalXOffset) {
			if (x_comp > totCellLength) {
				positions[i] = positions[i] - totCellLength ;
			}
			sendBuffers[Communicator::RIGHT].push_back(positions[i]) ;
			sendBuffers[Communicator::RIGHT].push_back(positions[i+1]) ;
			sendBuffers[Communicator::RIGHT].push_back(velocities[i]) ;
			sendBuffers[Communicator::RIGHT].push_back(velocities[i+1]) ;
			hasLeftProcess = true ;
		} else if (y_comp < bGlobalYOffset) {
			if (y_comp < 0) {
				positions[i+1] = totCellLength + positions[i+1] ;
			}
			sendBuffers[Communicator::BOT].push_back(positions[i]) ;
			sendBuffers[Communicator::BOT].push_back(positions[i+1]) ;
			sendBuffers[Communicator::BOT].push_back(velocities[i]) ;
			sendBuffers[Communicator::BOT].push_back(velocities[i+1]) ;
			hasLeftProcess = true ;
		} else if (y_comp > tGlobalYOffset) {
			std::cout << y_comp << std::endl;
			if (y_comp > totCellLength) {
				positions[i+1] = positions[i+1] - totCellLength ;
			}
			sendBuffers[Communicator::TOP].push_back(positions[i]) ;
			sendBuffers[Communicator::TOP].push_back(positions[i+1]) ;
			sendBuffers[Communicator::TOP].push_back(velocities[i]) ;
			sendBuffers[Communicator::TOP].push_back(velocities[i+1]) ;
			hasLeftProcess = true ;
		}
		if (hasLeftProcess) {
			std::iter_swap(positions.begin()+i+1,positions.end()-1);
  			positions.pop_back();
			std::iter_swap(positions.begin()+i,positions.end()-1);
  			positions.pop_back();
			std::iter_swap(velocities.begin()+i+1,velocities.end()-1);
  			velocities.pop_back();
			std::iter_swap(velocities.begin()+i,velocities.end()-1);
  			velocities.pop_back();
			forces.pop_back() ;
			forces.pop_back() ;
		} 
	}
	// Adds number of transitions to beginning of each send buffer. This is used to unpack later //
	for (int i = 0 ; i < 8 ; ++i) {
		sendBuffers[i].insert(sendBuffers[i].begin(),sendBuffers[i].size()/4) ;
	}
}

/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  fillDxZoneBuffers
 *  Description:  Detects if particles are near boundaries and fills the buffers,
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
void System<DataType,Integrator,Potential,Communicator>::communicateInteractionsTransitions() {
	Communicator::communicateBoundaries(sendBuffers,recvBuffers) ;
	// Unpack data. First handle transitions then interactions. //
	for (int i = 0 ; i < 8 ; ++i) {
		// Extract number of transitions. //
		unsigned int numTransitions = recvBuffers[i][0] ;
		unsigned int j = 1 ; 
		if (numTransitions !=  0) {
			for (j = 1 ; j < (numTransitions*4)+1 ; j+=4) {
				positions.push_back(recvBuffers[i][j]) ;
				positions.push_back(recvBuffers[i][j+1]) ;
				velocities.push_back(recvBuffers[i][j+2]) ;
				velocities.push_back(recvBuffers[i][j+3]) ;
				forces.push_back(0) ;
				forces.push_back(0) ;
			}
		}
		// Extract interactions. //
		for (j = j ; j < recvBuffers[i].size() ; j+=2) {
			boundaryGhostPositions.push_back(recvBuffers[i][j]) ;
			boundaryGhostPositions.push_back(recvBuffers[i][j+1]) ;
		}
	}
}

/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  printSystem
 *    Arguments:  std::ostream & out - Output stream used for printing.
 *                int iter - Iteration number.
 *  Description:  Prints position and velocity information of particles to output stream.
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
void System<DataType,Integrator,Potential,Communicator>::printSystem(std::ostream & out, int iter) const {
	if (positions.size() != 0) {
		for (unsigned int i = 0 ; i < positions.size() ; i+=2) {
			out << positions[i] << " " << positions[i+1] << " " ;
			out << velocities[i] << " " << velocities[i+1] << " " ;
			out << iter << " " << totCellLength << std::endl ;
		}
	} else {
		out << "-- -- -- -- -- " << totCellLength << std::endl  ;
	}
}		/* -----  end of member function printSystem  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  printEnergies
 *    Arguments:  std::ostream & out - Output stream used for printing.
 *                int iter - Iteration number.
 *  Description:  Prints relative energy deviation.
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
void System<DataType,Integrator,Potential,Communicator>::printEnergies(std::ostream & out, int iter) const {
	out << (initialEnergy - totEnergy) << " " << iter  << std::endl ;
}		/* -----  end of member function printSystem  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  ~System
 *  Description:  Destroy the particles in the system.
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential, class Communicator> 
System<DataType,Integrator,Potential,Communicator>::~System() {
	;
}		/* -----  end of member function ~System  ----- */

#endif /* end of include guard: SYSTEM_HPP_PSH64OVQ */
