CC = gcc
RM = rm -rf
USR_DIR = /home/9yelin9/usr
OMP_DIR = /opt/mpi/gcc-4.8.5/openmpi-4.1.0/
#CFLAGS = -O2 -Wall -mcmodel=medium -I$(USR_DIR)/include/lapacke -I$(OMP_DIR)/include -lm -llapack -fopenmp
CFLAGS = -g -O2 -Wall -mcmodel=medium -I$(USR_DIR)/include/lapacke -I$(OMP_DIR)/include -lm -llapack -fopenmp
LDFLAGS = -mcmodel=medium -L$(USR_DIR)/lib64 -L$(OMP_DIR)/lib
LINKS = -lm -llapack -fopenmp
OBJS = ../hf3.o cao.o 
TARGET = cao

.PHONY: all clean dep

all : $(TARGET)

clean :
	$(RM) *.o
	$(RM) $(TARGET)

dep :
	$(CC) $(CFLAGS) -M $(OBJS:.o=.c)

$(TARGET) : $(OBJS)
	$(CC) $(LDLIBS) $(LDFLAGS) -o $@ $(OBJS) $(LINKS)

.SUFFIXES : .c .o

.c .o :
	$(CC) $(CFLAGS) -c $<
