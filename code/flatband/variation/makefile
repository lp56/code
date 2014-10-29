CC = gcc
CFLAGS = -Wall -std=c99 -O3
LDLIBS = -lm -lgsl -lgslcblas

all: main
	./main

main : main.o var.o

main.o : main.c 

clean : 
	rm -f main.o var.o main
