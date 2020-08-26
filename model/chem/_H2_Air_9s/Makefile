include ../../../.makefile-header

LFLAGSLOCAL =  
CCFLAGSLOCAL = $(CFLAGSCFD) -I../../../lib/include -I../../../.

####### Files

HEADERS =   chem.h jachimowski1988old.h jachimowski1988.h
SOURCES =   chem.c jachimowski1988old.c jachimowski1988.c
OBJECTS =   chem.o jachimowski1988old.o jachimowski1988.o
TARGETS  = ../_chem.a


####### Build rules

all: src 

src: $(OBJECTS) 

arc: 
	( rm -f $(TARGETS) ; ar cqs $(TARGETS) *.o )                        

cleanall:
	( cd doc ; make clean )
	( make clean )

	

# DO NOT DELETE THIS LINE -- make depend depends on it.
