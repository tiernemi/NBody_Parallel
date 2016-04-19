/*
 * =====================================================================================
 *
 *       Filename:  leapfrog.cpp
 *
 *    Description:  Source file for leapfrog integration scheme. Non templated functions.
 *
 *        Version:  1.0
 *        Created:  04/13/2016 12:55:18 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

#include "leapfrog.hpp"

// Default valyue for deltaT
double Leapfrog::deltaT = 0.05 ;

/* 
 * ===  MEMBER FUNCTION CLASS : Leapfrog  ==============================================
 *         Name:  setDeltaT
 *    Arguments:  double deltaTarg - The time interval for the integrator.
 *  Description:  Setter fucntion for static time interval for the leapfrog integration
 *                policy.
 * =====================================================================================
 */

void Leapfrog::setDeltaT(double deltaTarg) {
	Leapfrog::deltaT = deltaTarg ;
}		/* -----  end of member function setDeltaT  ----- */


