CC = gcc
RM = rm -rf
USR_DIR = /home/9yelin9/usr
OMP_DIR = /opt/mpi/gcc-4.8.5/openmpi-4.1.0/
GSL_DIR = /opt/gsl/
CFLAGS = -g -O2 -Wall -mcmodel=medium -I../include -I$(USR_DIR)/include/lapacke -I$(OMP_DIR)/include -I$(GSL_DIR)/include -fopenmp
LDFLAGS = -mcmodel=medium -L../lib -L$(USR_DIR)/lib64 -L$(OMP_DIR)/lib -L$(GSL_DIR)/lib -fopenmp 
LINKS = -lhf3 -lm -llapack -lgsl -lgslcblas
OBJS = hf3.o
TARGET = hf3

.PHONY: clean dep
.SUFFIXES : .c .o

$(TARGET) : $(OBJS)
	$(CC) $(LDLIBS) $(LDFLAGS) -o $@ $(OBJS) $(LINKS)

.c .o :
	$(CC) $(CFLAGS) -c $<

clean :
	$(RM) *.o
	$(RM) $(TARGET)

dep :
	$(CC) $(CFLAGS) -M $(OBJS:.o=.c) 
