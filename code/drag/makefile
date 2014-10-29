CC = gcc
CFLAGS = -Wall -std=c99 -O3
LDLIBS = -lm -lgsl -lgslcblas

all: main
	./main.o >state.dat

main.o : main.c

main.c : ham.o

ham.o : ham.c
