MAKE = make
DIRS = lib hf3

.PHONY: all clean dep
.SUFFIXES : .c .o

all :
	@for d in $(DIRS); \
	do \
		$(MAKE) -C $$d; \
	done

clean :
	@for d in $(DIRS); \
	do \
		$(MAKE) -C $$d clean; \
	done
