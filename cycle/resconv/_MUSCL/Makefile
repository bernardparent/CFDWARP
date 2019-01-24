include ../../../.makefile-header

LFLAGSLOCAL =  
CCFLAGSLOCAL = $(CFLAGSCFD) -I../../../lib/include -I../../../.

####### Files

HEADERS =   resconv.h 
SOURCES =   resconv.c 
OBJECTS =   resconv.o  
TARGETS  =  ../_resconv.a


####### Build rules

all: src 

src: $(OBJECTS) 

arc: 
	( rm -f $(TARGETS) ; ar cqs $(TARGETS) *.o )

cleanall: clean
	( cd doc ; make clean )


# DO NOT DELETE THIS LINE -- make depend depends on it.
