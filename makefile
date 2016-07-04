CC=mpicxx
CXXFLAGS= -Wall -O3 -std=c++11
DEPS= system.hpp leapfrog.hpp leonard_jones.hpp mpi_utils.hpp torus_communicator.hpp time_utils.hpp
OBJ = main.o leapfrog.o mpi_utils.o torus_communicator.o time_utils.o
LDFlags= -lm
SRCDIR = ./src/
INCDIR = ./inc/
BIN = ./bin/
TARGET = nbod

CXXFILES= $(shell ls $(SRCDIR)*.cpp | xargs -n1 basename)
CXXOBJS= $(CXXFILES:cpp=o)
CXXSRC= $(addprefix $(SRCDIR), $CXXFILES)
CXXOBJECTS=$(addprefix $(BIN), $(CXXOBJS))


all: $(BIN) $(CXXOBJECTS) 
	$(CC) $(CXXOBJECTS) -o $(TARGET) -I$(INCDIR)

$(BIN):
	mkdir $(BIN)

clean:
	rm -rf $(BIN) $(TARGET) out.*

animation:
	mpirun -n 6 ./nbod -i 5000 -n 40 -f "out.txt" -a
	gnuplot -e "load \"animation.gnu\"" --persist
	rm -f "out.txt"

.SECONDEXPANSION:
$(CXXOBJECTS): %.o: $$(addprefix $(SRCDIR), $$(notdir %)).cpp
	$(CC) -c $< $(CXXFLAGS) -I$(INCDIR) -o $@

