include ../../../.makefile-header

CCFLAGSLOCAL = -I../../include
LFLAGSLOCAL = -L../../ -lgridg -lsoap -lexm -lm

####### Files

SOURCES =   tecgrid2D.c tecgrid3D.c
OBJECTS =   tecgrid2D.o tecgrid3D.o
TARGETS =   tecgrid2D tecgrid3D 

####### Build rules


all: $(TARGETS)

tecgrid2D: tecgrid2D.o
	$(L) tecgrid2D.o $(LFLAGS) $(LFLAGSLOCAL) -o tecgrid2D

tecgrid3D: tecgrid3D.o
	$(L) tecgrid3D.o $(LFLAGS) $(LFLAGSLOCAL) -o tecgrid3D

# DO NOT DELETE THIS LINE -- make depend depends on it.
