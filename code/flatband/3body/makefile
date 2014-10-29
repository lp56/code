CC = gcc
CFLAGS = -Wall -std=c99 -O3
LDLIBS = -lm -lgsl -lgslcblas

all : eigs.dat

eigs.dat : 3body
	./3body > eigs.dat

3body : 3body.o

3body.o : 3body.c

clean :
	rm 3body.o 3body eigs.dat 
