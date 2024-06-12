include ../../../.makefile-header

LFLAGSLOCAL =  
CCFLAGSLOCAL = $(CFLAGSCFD) -I../../../lib/include -I../../../.

####### Files

HEADERS =   chem.h macheret2007.h parent2023b.h parent2023.h rodriguez2024.h
SOURCES =   chem.c macheret2007.c parent2023b.c parent2023.c rodriguez2024.c
OBJECTS =   chem.o macheret2007.o parent2023b.o parent2023.o rodriguez2024.o
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
