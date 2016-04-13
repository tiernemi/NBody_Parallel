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

template <class DataType, class Integrator, class Potential> 
class System {
 public:
	System() {} ;
	virtual ~System() ;
	void addParticle(const std::pair<DataType,DataType> &, const std::pair<DataType,DataType> &) ;
	void simulate(int) ;
	void simulate(int, std::ostream &) ;
	int numParticles() const ;
	void printSystem(std::ostream &, int iter) const ;
 private:
	void updatePositions() ;
	void updateVelocities() ;
	void updateForces() ;
 private:
	std::vector<DataType> positions ;
	std::vector<DataType> velocities ;
	std::vector<DataType> forces ;
} ;		/* -----  end of class System  ----- */

// TEMPLATED MEMBER FUNCTIONS //

/* 
 * ===  MEMBER FUNCTION CLASS : System  =================================================
 *         Name:  System
 *    Arguments:  const std::pair<Datatype,Datatype> & position - Position of new particle.
 *                const std::pair<Datatype,Datatype> & velocity - Velocity of new particle.
 *  Description:  Add a new particle to the system.
 * ======================================================================================
 */

template <class DataType, class Integrator, class Potential> 
void System<DataType,Integrator,Potential>::addParticle(const std::pair<DataType,DataType> & position, const std::pair<DataType,DataType> & velocity) {
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

template <class DataType, class Integrator, class Potential> 
void System<DataType,Integrator,Potential>::simulate(int numIters) {
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
	}
}		/* -----  end of member function simulate  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  simulate
 *    Arguments:  int numIters - Number of iterations to run simulation.
 *                std::ostream & ouput - Output stream to print to.
 *  Description:  Simulate the system for a set number of iterations.
 * ====================================================================================
 */

template <class DataType, class Integrator, class Potential> 
void System<DataType,Integrator,Potential>::simulate(int numIters, std::ostream & output) {
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
		printSystem(output,i) ;
		output << std::endl << std::endl ;
	}
}		/* -----  end of member function simulate  ----- */



/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  numParticles
 *      Returns:  The number of particles in the system.
 *  Description:  Getter function for the number of particles in the system.
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential> 
int System<DataType,Integrator,Potential>::numParticles() const {
	return (positions.size()/2) ;
}		/* -----  end of member function num  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  updatePositions
 *  Description:  Update the positions of the particles in the system.
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential> 
void System<DataType,Integrator,Potential>::updatePositions() {
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

template <class DataType, class Integrator, class Potential> 
void System<DataType,Integrator,Potential>::updateVelocities() {	
	for (unsigned int i = 0 ; i < positions.size() ; i+=2) {
		Integrator::updateVelocity(&velocities[i],&forces[i]) ;
	}
}		/* -----  end of member function updateVelocities  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ======================================
 *         Name:  updateForces
 *  Description:  Update the forces of the particles in the system.
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential> 
void System<DataType,Integrator,Potential>::updateForces() {	
	unsigned int numParts = numParticles() ;
	for (auto i = forces.begin(); i != forces.end(); ++i) {
		*i = 0 ;
	}
	for (unsigned int i = 0 ; i < numParts ; ++i) {
		for (unsigned int j = i+1 ; j < numParts ; ++j) {
			DataType x_comp = positions[2*j] - positions[2*i] ;
			DataType y_comp = positions[2*j+1] - positions[2*i+1] ;
			DataType distSqr = x_comp*x_comp + y_comp*y_comp ;
			DataType forceMag = Potential::calcForce(distSqr) ;
			DataType xForce = forceMag*x_comp ;
			DataType yForce = forceMag*y_comp ;
			forces[2*i] += xForce ;
			forces[2*i+1] += yForce ;
			forces[2*j] -= xForce ;
			forces[2*j+1] -= yForce ;
		}
	}
}		/* -----  end of member function updateForces  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  printSystem
 *    Arguments:  std::ostream & out - Output stream used for printing.
 *                int iter - Iteration number.
 *  Description:  Prints position and velocity information of particles to output stream.
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential> 
void System<DataType,Integrator,Potential>::printSystem(std::ostream & out, int iter) const {
	for (unsigned int i = 0 ; i < positions.size() ; i+=2) {
		out << positions[i] << " " << positions[i+1] << " " ;
		out << velocities[i] << " " << velocities[i+1] << " " ;
		out << iter << std::endl ;
	}
}		/* -----  end of member function printSystem  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : System  ================================================
 *         Name:  ~System
 *  Description:  Destroy the particles in the system.
 * =====================================================================================
 */

template <class DataType, class Integrator, class Potential> 
System<DataType,Integrator,Potential>::~System() {
	;
}		/* -----  end of member function ~System  ----- */

#endif /* end of include guard: SYSTEM_HPP_PSH64OVQ */
