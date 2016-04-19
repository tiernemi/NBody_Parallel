CC=g++
CFLAGS= -Wall -Og -std=c++11
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
	./nbod -n 100000 > "out.txt"
	gnuplot -e "load \"animation.gnu\"" --persist
	rm -f "out.txt"


