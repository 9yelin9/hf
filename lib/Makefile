CC = gcc
AR = ar rv
RM = rm -rf

USR_DIR = /home/9yelin9/.local

CFLAGS = -g -O2 -Wall -mcmodel=medium -I../include -I$(USR_DIR)/include -fopenmp
OBJS = tb.o hf.o com.o
LIBS = libhf.a

.PHONY: all clean dep
.SUFFIXES : .c .o

.c .o :
	$(CC) $(CFLAGS) -c $<

all : $(LIBS)
clean :
	$(RM) *.o
dep :
	$(CC) $(CFLAGS) -M $(OBJS:.o=.c) 

$(LIBS) : $(OBJS)
	$(AR) $@ $(OBJS)
