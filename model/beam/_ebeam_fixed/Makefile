include ../../../.makefile-header

LFLAGSLOCAL =  
CCFLAGSLOCAL = $(CFLAGSCFD) -I../../../lib/include -I../../../.

####### Files

HEADERS =   beam.h 
SOURCES =   beam.c 
OBJECTS =   beam.o 
TARGETS  = ../_beam.a


####### Build rules

all: src 

src: $(OBJECTS) 

arc: 
	( rm -f $(TARGETS) ; ar cqs $(TARGETS) *.o )                        

cleanall: clean
	

# DO NOT DELETE THIS LINE -- make depend depends on it.
