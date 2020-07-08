SHELL = /bin/sh
CFLAGS = -c
CC = cc -Wall
L = cc
LFLAGS = -lm


####### Files

HEADERS =
SOURCES =   proprietary.c
OBJECTS =   proprietary.o
TARGETS =   proprietary

####### Build rules

default:
	make all

all: src

src:
	$(CC) $(CFLAGS) proprietary.c
	$(L) proprietary.o $(LFLAGS) -o proprietary

clean:
	rm -f *.o core $(TARGETS)

cleanall:
	rm -f *.o core $(TARGETS)


# DO NOT DELETE THIS LINE -- make depend depends on it.
