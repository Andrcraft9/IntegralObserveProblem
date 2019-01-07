CC = g++
CFLAGS = -O3

all:
	$(CC) $(CFLAGS) -o prog main.cpp solver.cpp

