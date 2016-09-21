# Makefile

PROGRAM = funecm
OBJS = point.o double_add.o normal_add.o ecm.o scalar.o main.o
CC = icc
# CFLAGS = -mmic -openmp -O2
CFLAGS = -openmp -O2
LIBDIR = -L/home/project8/gmp/lib
LIBRARY = -lgmp -lrt

# gccによるコンパイル
# CC = gcc
# CFLAGS = -fopenmp -O2
# LIBRARY = -lgmp -lrt -lm

$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) -o $(PROGRAM) $^ $(LIBRARY) $(LIBDIR)

main.o: main.c
	$(CC) $(CFLAGS) -c $<
point.o: point.c
	$(CC) $(CFLAGS) -c $<
double_add.o: double_add.c
	$(CC) $(CFLAGS) -c $<
normal_add.o: normal_add.c
	$(CC) $(CFLAGS) -c $<
ecm.o: ecm.c
	$(CC) $(CFLAGS) -c $<
scalar.o: scalar.c
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	$(RM) $(PROGRAM) $(OBJS)

