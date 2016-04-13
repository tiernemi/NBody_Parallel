/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  Main function for 2-body LJ simulation.
 *
 *        Version:  1.0
 *        Created:  04/12/2016 11:21:31 AM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

#include <iostream>
#include "getopt.h"
#include "system.hpp"
#include "leapfrog.hpp"
#include "leonard_jones.hpp"
#include <cstdlib>
#include <tuple>

int main(int argc, char *argv[]) {

	// Get options. //
	int numIters = 1E3 ;
	double deltaT = 0.05 ;
	int choice ;
	while (1) {
		choice = getopt(argc, argv, "n:t:");	
		if (choice == -1)
			break;
		switch( choice ) {
			case 'n':
				numIters = std::atoi(optarg) ;
				break;
			case 't':
				deltaT = std::atof(optarg) ;
				break;
			default:
				/* Not sure how to get here... */
				return EXIT_FAILURE;
		}
	}
	
	// Set integration time step. //
	Leapfrog::setDeltaT(deltaT) ;
	// Create a system of particles integrated using leapfrog and subject to a LJ potential.
	System<double,Leapfrog,LeonardJones> nbodySystem ;
	// Add particles to the system. //
	nbodySystem.addParticle(std::pair<double,double>(1,1), std::pair<double,double>(0,0)) ;
	nbodySystem.addParticle(std::pair<double,double>(-1,-1), std::pair<double,double>(0,0)) ;
	// Simulate the system and print to terminal. //
	nbodySystem.simulate(numIters, std::cout) ;

	return EXIT_SUCCESS ;
}
