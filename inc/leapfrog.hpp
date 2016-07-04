#ifndef LEAPFROG_HPP_ESX4AUNN
#define LEAPFROG_HPP_ESX4AUNN

/*
 * =====================================================================================
 *
 *       Filename:  leapfrog.hpp
 *
 *    Description:  Leapfrog integration policy. 
 *
 *        Version:  1.0
 *        Created:  04/12/2016 03:41:21 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

/* 
 * ===  CLASS  =========================================================================
 *         Name:  Leapfrog
 *       Fields:  double deltaT - Time interval used to update position and velocity.
 *  Description:  Leapfrog integration policy. This integration scheme uses a verlet-esque
 *                step to update positions and velocities. First the velocity is stepped
 *                back and then the x and v are a half step out of sync for each update.
 *                This is time reversible and as a result conserves energy.
 * =====================================================================================
 */

class Leapfrog {
 public:
	static void setDeltaT(double) ;
	template <class DataType>
	static inline void initialise(DataType *, DataType *, DataType *) ;
	template <class DataType>
	static inline void updatePosition(DataType *, DataType *) ;
	template <class DataType>
	static inline void updateVelocity(DataType *, DataType *) ;
 private:
	static double deltaT ;
} ;		/* -----  end of class Leapfrog  ----- */

// TEMPLATED MEMBER FUNCTIONS //

/* 
 * ===  MEMBER FUNCTION CLASS : Leapfrog  ==============================================
 *         Name:  initialise
 *    Arguments:  Datatype * position - A pointer to a 2D position array.
 *                Datatype * velocity - A pointer to a 2D velocity array.
 *                Datatype * force - A pointer to a 2D force array.
 *  Description:  Sets up the integrator by pre-conditioning positions,velocities and forces.
 * =====================================================================================
 */

template <class DataType>
inline void Leapfrog::initialise(DataType * position, DataType * velocity, DataType * force) {
	velocity[0] -= force[0] * deltaT/2.0 ;
	velocity[1] -= force[1] * deltaT/2.0 ;
}		/* -----  end of member function initialise  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : Leapfrog  =============================================
 *         Name:  updatePosition
 *    Arguments:  Datatype * position - A pointer to a 2D position array.
 *                Datatype * velocity - A pointer to a 2D velocity array.
 *  Description:  Updates the position using the leapfrog scheme.
 * =====================================================================================
 */

template <class DataType>
inline void Leapfrog::updatePosition(DataType * position, DataType * velocity) {	
	position[0] += velocity[0] * deltaT ;
	position[1] += velocity[1] * deltaT ;
}		/* -----  end of member function initialise  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : Leapfrog  =============================================
 *         Name:  updateVelocity
 *    Arguments:  Datatype * velocity - A pointer to a 2D velocity array.
 *                Datatype * force - A pointer to a 2D force array.
 *  Description:  Updates the velocity using the leapfrog scheme.
 * =====================================================================================
 */

template <class DataType>
inline void Leapfrog::updateVelocity(DataType * velocity, DataType * force) {
	velocity[0] += force[0] * deltaT ;
	velocity[1] += force[1] * deltaT ;
}		/* -----  end of member function initialise  ----- */

#endif /* end of include guard: LEAPFROG_HPP_ESX4AUNN */
