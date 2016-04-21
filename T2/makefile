CC=g++
CFLAGS= -Wall -O3 -std=c++11 -lrt -lm
DEPS= system.hpp leapfrog.hpp leonard_jones.hpp time_utils.hpp
OBJ = main.o leapfrog.o time_utils.o
LDFlags= -lm -lrt

nbod: $(OBJ)
		$(CC) -o $@ $^ $(CFLAGS)

leapfrog.o: leapfrog.cpp leapfrog.hpp
		$(CC) -c -o $@ $< $(CFLAGS)

time_utils.o: time_utils.cpp time_utils.hpp
		$(CC) -c -o $@ $< $(CFLAGS)

main.o: main.cpp $(DEPS)
		$(CC) -c -o $@ $< $(CFLAGS)


clean:
	rm -f $(OBJ) nbod

animation:
	./nbod -i 40000 -n 10 -f "out.txt" -a
	gnuplot -e "load \"animation.gnu\"" --persist
	rm -f "out.txt"

plot:
	./nbod -i 6000 -n 40 -e
	gnuplot -e "load \"plot.gnu\"" --persist
	rm -f "energy.txt"


