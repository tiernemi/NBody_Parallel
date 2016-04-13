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
} ;		/* -----  end of class LeonardJones  ----- */

// TEMPLATED MEMBER FUNCTIONS //

/* 
 * ===  MEMBER FUNCTION CLASS : LeonardJones  =========================================
 *         Name:  function
 *    Arguments:  const Datatype & distanceSqr - The distance between two particles squared.
 *      Returns:  The magnitude of the LJ force.
 *  Description:  The force calculation function for this particular type of potential.
 * =====================================================================================
 */

template <class DataType>
inline DataType LeonardJones::calcForce(const DataType & distanceSqr) {
	DataType r8 = std::pow(distanceSqr,-4) ;
	return r8 - std::pow(r8,3) ;
}		/* -----  end of member function function  ----- */

#endif /* end of include guard: POTENTIAL_HPP_36LKUZMX */
