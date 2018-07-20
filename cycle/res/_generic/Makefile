include ../../../.makefile-header

LFLAGSLOCAL =  
CCFLAGSLOCAL = $(CFLAGSCFD) -I../../../lib/include -I../../../.

####### Files

HEADERS =   res.h 
SOURCES =   res.c 
OBJECTS =   res.o  
TARGETS  =  ../_res.a


####### Build rules

all: src 

src: $(OBJECTS) 

arc: 
	( rm -f $(TARGETS) ; ar cqs $(TARGETS) *.o )

cleanall: clean

# DO NOT DELETE THIS LINE -- make depend depends on it.
