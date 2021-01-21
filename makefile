CC = g++
CFLAGS = -std=c++14 -Wall -Wextra -Werror -O3
EXEC = solver.out

all: clean
all: $(EXEC)

debug: CFLAGS += -DDEBUG
debug: clean
debug: $(EXEC)

$(EXEC): maphMain.o maphSat.o
	$(CC) $(CFLAGS) -o maph.out maphMain.o maphSat.o

maphMain.o: maphMain.cpp maphSat.hpp
	$(CC) $(CFLAGS) -c maphMain.cpp

maphSat.o: maphSat.cpp maphSat.hpp
	$(CC) $(CFLAGS) -c maphSat.cpp

clean:
	rm -f maph.out *.o
