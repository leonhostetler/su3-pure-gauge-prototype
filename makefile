#-------------------------------------------------#
#            Makefile for UNIX Systems            #
#             using a GNU C compiler              #
#-------------------------------------------------#
# The Makefile is a set of instructions for the
# make utility. It tells make how to compile a
# modular program using gcc.

# The following two lines are macros, e.g., every
# occurrence of $(CC) gets replaced by gcc
CC=gcc
CFLAGS=-g -Wall

# Description of compiler flags
# -g	adds debugging information to the executable
# -Wall	turns on most compiler warnings

# The target for the make utility
all: simulation

# Combine the object files using gcc into a program called "simulation"
# The "-lm" must be added when linking the math library
simulation: main.o su3matrices.o complex.o
	$(CC) $(CFLAGS) -o simulation main.o su3matrices.o complex.o -lm

# Make the object files by combining c and header files
main.o: main.c su3matrices.h
su3matrices.o: su3matrices.c su3matrices.h
complex.o: complex.c complex.h

# The command that is run if "make clean" is invoked
clean:
	rm -f simulation main.o su3matrices.o complex.o
