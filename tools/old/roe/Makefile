include ../../.makefile-header

CCFLAGSLOCAL = -I../../lib/include
LFLAGSLOCAL = -L../../lib -lm -lextmaths -lextstdio

####### Files

HEADERS =       
SOURCES =	main.c   
OBJECTS =	main.o     
TARGETS  =	roe

####### Build rules

all: roe

roe: $(OBJECTS) 
	$(L) $(OBJECTS) -o $(TARGETS) $(LFLAGS) $(LFLAGSLOCAL)

# DO NOT DELETE THIS LINE -- make depend depends on it.
