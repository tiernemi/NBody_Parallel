//..................................................................................//
______               _ _      _        _   _       ______           _       
| ___ \             | | |    | |      | \ | |      | ___ \         | |      
| |_/ /_ _ _ __ __ _| | | ___| |      |  \| |______| |_/ / ___   __| |_   _ 
|  __/ _` | '__/ _` | | |/ _ \ |      | . ` |______| ___ \/ _ \ / _` | | | |
| | | (_| | | | (_| | | |  __/ |      | |\  |      | |_/ / (_) | (_| | |_| |
\_|  \__,_|_|  \__,_|_|_|\___|_|      \_| \_/      \____/ \___/ \__,_|\__, |
                                                                       __/ |
                                                                      |___/ 
//.................................................................................//

The following code contains a generic implementation of a parallel N-body solver. 
Integration methods, potentials and the communication regime are specified at compile time using
policy templates. The current implementation uses a leapfrog integrator, LJ potential and a torus
communicator (Periodic boundaries). Interprocess communication is done using MPI. Animation
is done using a simple gnuplot script. Program will output final positions of particles after
simulation or can generate an animation of the simulation using the -a flag in conjunction with the
animation.gnu gnuplot script.

REQUIREMENTS :
	OPENMPI 1.8.6 or greater : https://www.open-mpi.org/software/ompi/v1.10/
	GNUPLOT 5.0.3 or greater
        GCC 4.9.3 or greater with MPICXX.
	Some way to run MPI such as mpirun or mpiexec.

COMPILATION :
	Compilation just requires running "make" in the directory containing the makefile.

EXECUTION :
	Use mpirun -n [NPROC] ./nbod [OPTIONS]
	-i : Number of iterations [DEFAULT = 1E3] 
	-n : Number of particles [DEFAULT = 10] 
	-l : Simulation cell length [DEFAULT = 10] 
	-t : Width of timestep [DEFAULT = 0.001] 
	-f : Output file name [DEFAULT = DISABLED] 
	-a : Animation mode [DEFAULT = false] 
	-b : Benchmarking mode [DEFAULT = false] 
	-h : Help Message 
        ./nbod -h displays help information.
	
ANIMATION :
	Use make animation to run a trial animation specified in the makefile. Otherwise
        run mpirun -n [NPROC] ./nbod -a -f "data/out.txt" [OPTIONS] and then run "load animation.gnu"
	inside gnuplot. It's easier just to edit the parameters in the makefile.
