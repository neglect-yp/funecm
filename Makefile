# Makefile

PROGRAM = funecm
OBJS = atkin.o point.o double_add.o normal_add.o ecm.o scalar.o main.o
CC = icc
LIBRARY = -lgmp -lrt

ifeq ($(ARCHITECTURE),Xeonphi) # for Xeonphi
	CFLAGS = -mmic -openmp -O2
	LIBDIR = -L/usr/local/lib
else # for Xeon
	CFLAGS = -openmp -O2
	LIBDIR = -L/home/project8/gmp/lib
endif

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
atkin.o: atkin.c
	$(CC) $(CFLAGS) -c $<

all: clean $(PROGRAM)

clean:
	-rm -f $(PROGRAM) $(OBJS)

