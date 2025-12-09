CC = gcc
CFLAGS = -Wall -g -O3
ifeq ($(FULL),1)
CPPFLAGS += -DFULL_MULTIPLICITY
endif

TARGETS = keyrec
OBJS = main.o poly.o spectrum.o

all: $(TARGETS)

keyrec: $(OBJS)
	$(CC) $(LDFLAGS) $^ -o $@

clean:
	/bin/rm $(OBJS)

veryclean: clean
	/bin/rm $(TARGETS)

main.o: poly.h spectrum.h
poly.o: poly.h
