include ../../../.makefile-header

LFLAGSLOCAL =  
CCFLAGSLOCAL = $(CFLAGSCFD) -I../../../lib/include -I../../../.

####### Files

HEADERS =   thermo.h enthalpy.h electronimpact.h
SOURCES =   thermo.c enthalpy.c electronimpact.c
OBJECTS =   thermo.o  enthalpy.o  electronimpact.o
TARGETS  = ../_thermo.a


####### Build rules

all: src 

src: $(OBJECTS) 

arc: 
	( rm -f $(TARGETS) ; ar cqs $(TARGETS) *.o )                        

cleanall: clean

# DO NOT DELETE THIS LINE -- make depend depends on it.
