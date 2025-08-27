include ../../../.makefile-header

LFLAGSLOCAL =  
CCFLAGSLOCAL = $(CFLAGSCFD) -I../../../lib/include -I../../../.

####### Files

HEADERS =   ts.h 
SOURCES =   ts.c 
OBJECTS =   ts.o  
TARGETS  =  ../_ts.a


####### Build rules

all: src 

src: $(OBJECTS) 

arc: 
	( rm -f $(TARGETS) ; ar cqs $(TARGETS) *.o )


# DO NOT DELETE THIS LINE -- make depend depends on it.
