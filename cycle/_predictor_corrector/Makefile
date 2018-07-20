include ../../.makefile-header

LFLAGSLOCAL =  
CCFLAGSLOCAL = $(CFLAGSCFD) -I../../lib/include -I../../.

####### Files

HEADERS =   cycle.h 
SOURCES =   cycle.c 
OBJECTS =   cycle.o  
TARGETS  =  ../_cycle.a


####### Build rules

all: src 

src: $(OBJECTS) 

arc: 
	ar -x ../ts/_ts.a  
	ar -x ../tsemf/_tsemf.a  
	ar -x ../res/_res.a  
	ar -x ../restime/_restime.a  
	ar -x ../resconv/_resconv.a  
	ar -x ../ressource/_ressource.a  
	ar -x ../share/cycle_share.a
	( rm -f $(TARGETS) ; ar cqs $(TARGETS) *.o )

# DO NOT DELETE THIS LINE -- make depend depends on it.
