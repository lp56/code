CC = gcc
CFLAGS = -Wall -std=c99 -O3
LDLIBS = -lm -lgsl -lgslcblas

all : eigs.dat

eigs.dat : 2body
	./2body > eigs.dat

2body : 2body.o

2body.o : 2body.c

clean :
	rm -f 2body.o 2body eigs.dat
