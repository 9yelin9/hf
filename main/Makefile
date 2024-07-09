CC = gcc
RM = rm -rf

USR_DIR = /home/9yelin9/.local
OPENBLAS_DIR = /opt/openblas/gcc-11.4.1/0.3.26

CFLAGS = -g -O2 -Wall -mcmodel=medium -I../include -I$(USR_DIR)/include -I$(OPENBLAS_DIR)/include -fopenmp
LDFLAGS = -L../lib -L$(USR_DIR)/lib -L$(OPENBLAS_DIR)/lib -fopenmp 
#LINKS = -lhf -lz -lm -lopenblas -lgsl -lgslcblas -lhdf5
LINKS = -lhf -lz -lm -llapack -lblas -lopenblas -lhdf5

OBJS = tb.o hf.o
TARGETS = ../tb ../hf

.PHONY : all clean dep
.SUFFIXES : .c .o

.c .o :
	$(CC) $(CFLAGS) -c $<

all : $(TARGETS)
clean :
	$(RM) *.o
	$(RM) $(TARGETS)
dep :
	$(CC) $(CFLAGS) -M $(OBJS:.o=.c) 

../tb : tb.o
	$(CC) $(LDLIBS) $(LDFLAGS) -o $@ tb.o $(LINKS)
../hf : hf.o
	$(CC) $(LDLIBS) $(LDFLAGS) -o $@ hf.o $(LINKS)
