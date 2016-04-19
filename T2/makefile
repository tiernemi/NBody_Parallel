CC=g++
CFLAGS= -Wall -O3 -std=c++11
DEPS= system.hpp leapfrog.hpp leonard_jones.hpp
OBJ = main.o leapfrog.o
LDFlags= -lm

nbod: $(OBJ)
		$(CC) -o $@ $^ $(CFLAGS)

leapfrog.o: leapfrog.cpp leapfrog.hpp
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
	./nbod -i 400000 -n 10 -e
	gnuplot -e "load \"plot.gnu\"" --persist
	rm -f "energy.txt"


