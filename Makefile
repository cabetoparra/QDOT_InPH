#PREDEFINED VARIABLES
CC = gcc
FF = gfortran

CFLAGS = -g -c -O3 -fopenmp -w -L/usr/local/lib 
LFLAGS = -lm -O3 -fopenmp -lgsl -lgslcblas -llapack

PROGRAM = prog-tester


MODULOSC = $(PROGRAM).c libraries.c
MODULOSO = $(PROGRAM).o libraries.o 

run:	compile
	@./$(PROGRAM) 

compile:	clean
	@$(CC) $(CFLAGS) $(MODULOSC)
	@$(CC) -o $(PROGRAM) $(MODULOSO) $(LFLAGS)

clean:
	@rm -rf *.o $(PROGRAM) *.c~ *.dat~ *~ *.log *.mod
