CC = gcc
CFLAGS = -Wall -std=c99 -O3
LDLIBS = -lm -lgsl -lgslcblas

all: main
	./main > dat.dat

main : main.o 2bhw.o alpha.o mtr.o

main.o : main.c 

clean : 
	rm -f main.o 2bhw.o alpha.o mtr.o main
