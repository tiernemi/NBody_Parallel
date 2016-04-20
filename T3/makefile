CC=mpicxx
CFLAGS= -Wall -O3 -std=c++11
DEPS= system.hpp leapfrog.hpp leonard_jones.hpp mpi_utils.hpp torus_communicator.hpp
OBJ = main.o leapfrog.o mpi_utils.o torus_communicator.cpp
LDFlags= -lm

nbod: $(OBJ)
		$(CC) -o $@ $^ $(CFLAGS)

leapfrog.o: leapfrog.cpp leapfrog.hpp
		$(CC) -c -o $@ $< $(CFLAGS)


mpi_utils.o: mpi_utils.cpp mpi_utils.hpp
		$(CC) -c -o $@ $< $(CFLAGS)
		
torus_communicator.o: torus_communicator.cpp torus_communicator.hpp
		$(CC) -c -o $@ $< $(CFLAGS)

main.o: main.cpp $(DEPS)
		$(CC) -c -o $@ $< $(CFLAGS)


clean:
	rm -f $(OBJ) nbod

animation:
	mpirun -n 4 ./nbod -i 2000 -n 10 -f "out.txt" -a -t 0.001
	gnuplot -e "load \"animation.gnu\"" --persist
	rm -f "out.txt"

plot:
	mpirun -n 4 ./nbod -i 400000 -n 10 -e
	gnuplot -e "load \"plot.gnu\"" --persist
	rm -f "energy.txt"


