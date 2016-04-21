/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  Main function for N-body LJ simulation. This system has periodic
 *                  boundaries and a uniformly distributed set of particles as an initial
 *                  state.
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
#include <fstream>
#include <tuple>
#include "time_utils.hpp"

int main(int argc, char *argv[]) {

	// Get options. //
	int numIters = 1E3 ; // Number of iterations to run for. 
	int numParticles = 10 ;
	double deltaT = 0.001 ; // The size of the timestep
	double cellLength = 10.0f ; // The length of the cell.
	bool trackEnergy = false ; // Flag for printing energy.
	bool animation = false ; // Flag for creating animation output.
	std::string filename = "" ; // Filename string.
	bool benching = false ;
	int choice ;
	while (1) {
		choice = getopt(argc, argv, "i:n:t:l:ef:ab");	
		if (choice == -1)
			break;
		switch( choice ) {
			case 'i':
				numIters = std::atoi(optarg) ;
				break;
			case 'n':
				numParticles = std::atoi(optarg) ;
				break;
			case 'l':
				cellLength = std::atof(optarg) ;
				break;
			case 't':
				deltaT = std::atof(optarg) ;
				break;
			case 'e':
				trackEnergy = true ;
				break;
			case 'f':
				filename = optarg ;
				break;
			case 'a':
				animation = true ;
				break;
			case 'b':
				benching  = true ;
				break;
			default:
				/* Not sure how to get here... */
				return EXIT_FAILURE;
		}
	}
	
	if (benching) {
		startClock() ;
	}
	// Set integration time step. //
	Leapfrog::setDeltaT(deltaT) ;
	// Create a system of particles integrated using leapfrog and subject to a LJ potential.
	System<double,Leapfrog,LeonardJones> nbodySystem(cellLength) ;
	// Add particles to the system using a uniform distribution for positions. //
	nbodySystem.randomInitialise(15071992, numParticles) ;
	
	if (animation) {
		// If no file name passed then print to cout. //
		if (filename.compare("") == 0) {
			if (trackEnergy) {
				std::ofstream en("energy.txt") ;
				nbodySystem.simulateEnergies(numIters, std::cout, en) ;
				en.close() ;
			} else {
				nbodySystem.simulate(numIters, std::cout) ;
			}
		} 
		// Else print to filename.
		else {
			std::ofstream posOutput(filename) ;
			if (trackEnergy) {
				std::ofstream en("energy.txt") ;
				nbodySystem.simulateEnergies(numIters, posOutput, en) ;
				en.close() ;
			} else {
				nbodySystem.simulate(numIters, posOutput) ;
			}
			posOutput.close() ;
		}
	} else {
	// If no file name passed then print to cout. //
		if (filename.compare("") == 0) {
			if (trackEnergy) {
				std::ofstream en("energy.txt") ;
				nbodySystem.simulateEnergies(numIters, en) ;
				nbodySystem.printSystem(std::cout,numIters) ;
				en.close() ;
			} else {
				nbodySystem.simulate(numIters) ;
				nbodySystem.printSystem(std::cout,numIters) ;
			}
		} 
		// Else print to filename.
		else {
			std::ofstream posOutput(filename) ;
			if (trackEnergy) {
				std::ofstream en("energy.txt") ;
				nbodySystem.simulateEnergies(numIters, en) ;
				nbodySystem.printSystem(posOutput,numIters) ;
				en.close() ;
			} else {
				nbodySystem.simulate(numIters, posOutput) ;
				nbodySystem.printSystem(posOutput,numIters) ;
			}
			posOutput.close() ;
		}

	}

	if (benching) {
		stopClock() ;
		std::cout << getElapsedTime() << std::endl ;
	}

	return EXIT_SUCCESS ;
}
