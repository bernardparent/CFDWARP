include ../../../.makefile-header

LFLAGSLOCAL =  
CCFLAGSLOCAL = $(CFLAGSCFD) -I../../../lib/include -I../../../.

####### Files

HEADERS =   chem.h dunnkang1973.h park1993.h boyd2007.h
SOURCES =   chem.c dunnkang1973.c park1993.c boyd2007.c
OBJECTS =   chem.o dunnkang1973.o park1993.o boyd2007.o
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
