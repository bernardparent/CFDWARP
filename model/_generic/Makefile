include ../../.makefile-header

LFLAGSLOCAL =  
CCFLAGSLOCAL = $(CFLAGSCFD) -I../../lib/include -I../../.

####### Files

HEADERS =   model.h 
SOURCES =   model.c
OBJECTS =   model.o
TARGETS  =  ../_model.a


####### Build rules

all: src 

src: $(OBJECTS) 

arc: 
	ar -x ../fluid/_fluid.a  
	ar -x ../chem/_chem.a  
	ar -x ../thermo/_thermo.a  
	ar -x ../metrics/_metrics.a
	ar -x ../emfield/_emfield.a
	ar -x ../beam/_beam.a
	ar -x ../share/model_share.a
	( rm -f $(TARGETS) ; ar cqs $(TARGETS) *.o )

cleanall: clean

# DO NOT DELETE THIS LINE -- make depend depends on it.
