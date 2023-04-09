CC = gcc
RM = rm -rf
USR_DIR = /home/9yelin9/.local
OMP_DIR = /opt/mpi/gcc-4.8.5/openmpi-4.1.0
GSL_DIR = /opt/gsl
HDF_DIR = /opt/hdf5/intel-2021.4.0/1.12.0
CFLAGS = -g -O2 -Wall -mcmodel=medium -I../include -I$(USR_DIR)/include -I$(OMP_DIR)/include -I$(GSL_DIR)/include -I$(HDF_DIR)/include -fopenmp
LDFLAGS = -L../lib -L$(USR_DIR)/lib -L$(OMP_DIR)/lib -L$(GSL_DIR)/lib -L$(HDF_DIR)/lib -fopenmp 
LINKS = -lhf3 -lz -lm -lopenblas -lgsl -lgslcblas -lhdf5
OBJS = hf3.o
TARGETS = hf3

.PHONY: all clean dep
.SUFFIXES : .c .o

.c .o :
	$(CC) $(CFLAGS) -c $<

all : $(TARGETS)

in : in.o
	$(CC) $(LDLIBS) $(LDFLAGS) -o $@ in.o $(LINKS)

hf3 : hf3.o
	$(CC) $(LDLIBS) $(LDFLAGS) -o $@ hf3.o $(LINKS)

out : out.o
	$(CC) $(LDLIBS) $(LDFLAGS) -o $@ out.o $(LINKS)

clean :
	$(RM) *.o
	$(RM) $(TARGET)

dep :
	$(CC) $(CFLAGS) -M $(OBJS:.o=.c) 
