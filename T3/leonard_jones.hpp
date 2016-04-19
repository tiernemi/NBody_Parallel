#ifndef POTENTIAL_HPP_36LKUZMX
#define POTENTIAL_HPP_36LKUZMX

/*
 * =====================================================================================
 *
 *       Filename:  leonard_jones.hpp
 *
 *    Description:  Leonard Jones Potential policy.
 *
 *        Version:  1.0
 *        Created:  04/12/2016 05:44:45 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

#include <cmath>

/* 
 * ===  CLASS  =========================================================================
 *         Name:  LeonardJones
 *  Description:  Leonard-Jones potential policy. This particular potential is of the
 *                form (r^2)^(-4) - (r^2)^(-7). A short ranges it's repulsive and at long
 *                ranges it's attractive. All constants set to 1.
 * =====================================================================================
 */

class LeonardJones {
 public:
	 template <class DataType>
	 static inline DataType calcForce(const DataType & distance) ;
	 template <class DataType>
	 static inline DataType calcEnergy(const DataType & distance) ;
 public:
	 static constexpr double cutOffDistSqr = 5 ;
} ;		/* -----  end of class LeonardJones  ----- */

// TEMPLATED MEMBER FUNCTIONS //

/* 
 * ===  MEMBER FUNCTION CLASS : LeonardJones  =========================================
 *         Name:  calcForce
 *    Arguments:  const Datatype & distanceSqr - The distance between two particles squared.
 *      Returns:  The magnitude of the LJ force.
 *  Description:  The force calculation function for this particular type of potential.
 * =====================================================================================
 */

template <class DataType>
inline DataType LeonardJones::calcForce(const DataType & distanceSqr) {
	return std::pow(distanceSqr,-4) - std::pow(distanceSqr,-7) ;
}		/* -----  end of member function function  ----- */

/* 
 * ===  MEMBER FUNCTION CLASS : LeonardJones  =========================================
 *         Name:  calcEnergy
 *    Arguments:  const Datatype & distanceSqr - The distance between two particles squared.
 *      Returns:  The magnitude of the LJ potential.
 *  Description:  The energy calculation function for this particular type of potential.
 * =====================================================================================
 */

template <class DataType>
inline DataType LeonardJones::calcEnergy(const DataType & distanceSqr) {
	return -(1/6.)*std::pow(distanceSqr,-3) + std::pow(distanceSqr,-6)/12. ;
}		/* -----  end of member function function  ----- */



#endif /* end of include guard: POTENTIAL_HPP_36LKUZMX */
