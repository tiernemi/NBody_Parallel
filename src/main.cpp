/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  Main function for N-body LJ simulation. This system has periodic
 *                  boundaries and a uniformly distributed set of particles as an initial
 *                  state. Program outputs position data to be used for plotting or
 *                  to generate animations.
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
#include <cstdlib>
#include <fstream>
#include "mpi.h"

#include "mpi_utils.hpp"
#include "torus_communicator.hpp"
#include "time_utils.hpp"
#include "system.hpp"
#include "leapfrog.hpp"
#include "leonard_jones.hpp"

static void displayHelpMessage() ;

int main(int argc, char *argv[]) {

	MPIUtils::initMPI(argc, argv) ;
	// Get options. //
	int numIters = 1000 ; // Number of iterations to run for. 
	int numParticles = 10 ;
	double deltaT = 0.001 ; // The size of the timestep
	double cellLength = 10.0f ; // The length of the cell.
	bool animation = false ; // Flag for creating animation output.
	bool benching = false ;
	std::string filename = "" ; // Filename string.
	int choice ;
	while (1) {
		choice = getopt(argc, argv, "i:n:t:l:ef:abh");	
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
			case 'f':
				filename = optarg ;
				break;
			case 'a':
				animation = true ;
				break;
			case 'b':
				benching = true ;
				break;
			case 'h':
				displayHelpMessage() ;
				return EXIT_SUCCESS ;
				break;
			default:
				return EXIT_FAILURE ;
		}
	}
	if (benching && MPIUtils::rank == 0) {
		startClock() ;
	}
	
	// Set integration time step. //
	Leapfrog::setDeltaT(deltaT) ;
	// Create a system of particles integrated using leapfrog and subject to a LJ potential.
	System<double,Leapfrog,LeonardJones,TorusCommunicator> nbodySystem(cellLength) ;
	// Add particles to the system using a uniform distribution for positions. //
	nbodySystem.randomInitialise(151992*MPIUtils::rank, numParticles) ;

	// If animation enabled, continuously stream position data to files/cout . //
	if (animation) {
		// If no file name passed then print to cout. //
		if (filename.compare("") == 0) {
			std::cout << "X-Co   | Y-Co   | X-Vel   | Y-Vel   | Iteration Number | BOXDIM" << std::endl ;
			nbodySystem.simulate(numIters, std::cout) ;
		} 
		// Else print to filename. //
		else {
			filename += std::to_string(MPIUtils::rank) ;
			std::ofstream posOutput(filename) ;
			nbodySystem.simulate(numIters, posOutput) ;
			posOutput.close() ;
		}
	} else {
	// If no file name passed then print to cout. //
		if (filename.compare("") == 0) {
			std::cout << "X-Co   | Y-Co   | X-Vel   | Y-Vel   | Iteration Number | BOXDIM" << std::endl ;
			nbodySystem.simulate(numIters) ;
			nbodySystem.printSystem(std::cout,numIters) ;
		} 
		// Else print to filename.
		else {
			filename += std::to_string(MPIUtils::rank) ;
			std::ofstream posOutput(filename) ;
			nbodySystem.simulate(numIters, posOutput) ;
			nbodySystem.printSystem(posOutput,numIters) ;
			posOutput.close() ;
		}

	}

	if (benching && MPIUtils::rank == 0) {
		stopClock() ;
		std::cout << getElapsedTime() << std::endl;
	}

	MPIUtils::finaliseMPI() ;
	return EXIT_SUCCESS ;
}

static void displayHelpMessage() {
	printf("This program simulates a generalised parallel N-body system given some potential, integrator and communicator.\n") ;
	printf("The current implementation uses a leapfrog integrator, torus commuincator and leonard jones potential.\n") ;
	printf("Usage : mpirun -n [NPROC] ./nbod [OPTIONS] \n") ;
	printf("Arguments : \n") ;
	printf("-i : Number of iterations [DEFAULT = 1E3] \n") ;
	printf("-n : Number of particles [DEFAULT = 10] \n") ;
	printf("-l : Simulation cell length [DEFAULT = 10] \n") ;
	printf("-t : Width of timestep [DEFAULT = 0.001] \n") ;
	printf("-f : Output file name [DEFAULT = DISABLED] \n") ;
	printf("-a : Animation mode [DEFAULT = false] \n") ;
	printf("-b : Benchmarking mode [DEFAULT = false] \n") ;
	printf("-h : Help Message \n") ;
}
