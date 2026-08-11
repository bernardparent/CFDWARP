include ../../../.makefile-header

CCFLAGSLOCAL = -I../../include
LFLAGSLOCAL = -L../../ -lgridg -lsoap -lexm -lm

####### Files

SOURCES =   printgrid1D.c tecgrid1D.c tecgrid2D.c tecgrid3D.c
OBJECTS =   printgrid1D.o tecgrid1D.o tecgrid2D.o tecgrid3D.o
TARGETS =   printgrid1D tecgrid1D tecgrid2D tecgrid3D 

####### Build rules


all: $(TARGETS)

printgrid1D: printgrid1D.o
	$(L) printgrid1D.o $(LFLAGS) $(LFLAGSLOCAL) -o printgrid1D

tecgrid1D: tecgrid1D.o
	$(L) tecgrid1D.o $(LFLAGS) $(LFLAGSLOCAL) -o tecgrid1D

tecgrid2D: tecgrid2D.o
	$(L) tecgrid2D.o $(LFLAGS) $(LFLAGSLOCAL) -o tecgrid2D

tecgrid3D: tecgrid3D.o
	$(L) tecgrid3D.o $(LFLAGS) $(LFLAGSLOCAL) -o tecgrid3D

# DO NOT DELETE THIS LINE -- make depend depends on it.
