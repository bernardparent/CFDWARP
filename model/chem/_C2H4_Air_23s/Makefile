include ../../../.makefile-header

LFLAGSLOCAL =  
CCFLAGSLOCAL = $(CFLAGSCFD) -I../../../lib/include -I../../../.

####### Files

HEADERS =   chem.h zettervall2017.h zettervall2017b.h konnovzettervall.h
SOURCES =   chem.c zettervall2017.c zettervall2017b.c konnovzettervall.c
OBJECTS =   chem.o zettervall2017.o zettervall2017b.o konnovzettervall.o
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
