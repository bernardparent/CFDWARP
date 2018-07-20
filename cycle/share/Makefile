include ../../.makefile-header

LFLAGSLOCAL =  
CCFLAGSLOCAL = $(CFLAGSCFD) -I../../lib/include -I../../.

####### Files

HEADERS =   cycle_share.h res_share.h ts_share.h tsemf_share.h 
SOURCES =   cycle_share.c res_share.c ts_share.c tsemf_share.c
OBJECTS =   cycle_share.o res_share.o ts_share.o tsemf_share.o
TARGETS  =  cycle_share.a


####### Build rules

all: src 

src: $(OBJECTS) 

arc: 
	( rm -f $(TARGETS) ; ar cqs $(TARGETS) *.o )

cleanall: clean

# DO NOT DELETE THIS LINE -- make depend depends on it.
