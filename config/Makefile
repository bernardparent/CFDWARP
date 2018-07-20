SHELL = /bin/sh
CFLAGS = -c
CC = cc -Wall
L = cc
LFLAGS = -lm


####### Files

HEADERS =
SOURCES =   config.c
OBJECTS =   config.o
TARGETS =   config

####### Build rules

default:
	make all

all: config

config:
	$(CC) $(CFLAGS) config.c
	$(L) config.o $(LFLAGS) -o config	

configx:
	$(CC) $(CFLAGS) -DEXPERIMENTAL config.c
	$(L) config.o $(LFLAGS) -o config

configp:
	$(CC) $(CFLAGS) -DPROPRIETARY config.c
	$(L) config.o $(LFLAGS) -o config

configpx:
	$(CC) $(CFLAGS) -DPROPRIETARY -DEXPERIMENTAL config.c
	$(L) config.o $(LFLAGS) -o config

clean:
	rm -f *.o core $(TARGETS)
	rm -f applyconfig.sh

cleanall: clean
	( cd _proprietary ; make clean ) 


# DO NOT DELETE THIS LINE -- make depend depends on it.    
