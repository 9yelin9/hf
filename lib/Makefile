CC = gcc
AR = ar rv
RM = rm -rf
USR_DIR = /home/9yelin9/usr
OMP_DIR = /opt/mpi/gcc-4.8.5/openmpi-4.1.0/
GSL_DIR = /opt/gsl/
CFLAGS = -g -O2 -Wall -mcmodel=medium -I../include -I$(USR_DIR)/include/lapacke -I$(OMP_DIR)/include -I$(GSL_DIR)/include -fopenmp
LDFLAGS = -mcmodel=medium -L$(USR_DIR)/lib64 -L$(OMP_DIR)/lib -L$(GSL_DIR)/lib -fopenmp 
LINKS = -lm -llapack -lgsl -lgslcblas
OBJS = func.o path.o
LIB = libhf3.a

.PHONY: clean
.SUFFIXES : .c .o

$(LIB) : $(OBJS)
	$(AR) $@ $(OBJS)

.c .o :
	$(CC) $(CFLAGS) -c $<

clean :
	$(RM) *.o
