CC = gcc
AR = ar rv
RM = rm -rf
USR_DIR = /home/9yelin9/.local
OMP_DIR = /opt/mpi/gcc-4.8.5/openmpi-4.1.0
GSL_DIR = /opt/gsl
HDF_DIR = /opt/hdf5/intel-2021.4.0/1.12.0
CFLAGS = -g -O2 -Wall -mcmodel=medium -I../include -I$(USR_DIR)/include -I$(OMP_DIR)/include -I$(GSL_DIR)/include -I$(HDF_DIR)/include -fopenmp
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
