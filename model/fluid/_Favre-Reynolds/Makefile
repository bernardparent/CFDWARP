include ../../../.makefile-header

LFLAGSLOCAL =  
CCFLAGSLOCAL = $(CFLAGSCFD) -I../../../lib/include -I../../../.

####### Files

HEADERS =   fluid.h fluid_post.h fluid_conv.h fluid_source.h fluid_init.h fluid_bdry.h 
SOURCES =   fluid.c fluid_post.c fluid_conv.c fluid_source.c fluid_init.c fluid_bdry.c 
OBJECTS =   fluid.o fluid_post.o fluid_conv.o fluid_source.o fluid_init.o fluid_bdry.o 
TARGETS  =  ../_fluid.a


####### Build rules

all: src 

src: $(OBJECTS) 

arc: 
	( rm -f $(TARGETS) ; ar cqs $(TARGETS) *.o )                        

cleanall: clean

	

# DO NOT DELETE THIS LINE -- make depend depends on it.
