CC = gcc
AR = ar rv
RM = rm -rf
USR_DIR = /home/9yelin9/usr
OMP_DIR = /opt/mpi/gcc-4.8.5/openmpi-4.1.0/
GSL_DIR = /opt/gsl/
CFLAGS = -g -O2 -Wall -mcmodel=medium -I../include -I$(USR_DIR)/include/lapacke -I$(OMP_DIR)/include -I$(GSL_DIR)/include -fopenmp
LDFLAGS = -mcmodel=medium -L$(USR_DIR)/lib64 -L$(OMP_DIR)/lib -L$(GSL_DIR)/lib -fopenmp 
LINKS = -lm -llapack -lgsl -lgslcblas
OBJS = func.o
LIBS = libhf3.a

.PHONY: all clean dep
.SUFFIXES : .c .o

.c .o :
	$(CC) $(CFLAGS) -c $<

all : $(LIBS)

$(LIBS) : $(OBJS)
	$(AR) $@ $(OBJS)

clean :
	$(RM) *.o

dep :
	$(CC) $(CFLAGS) -M $(OBJS:.o=.c) 
