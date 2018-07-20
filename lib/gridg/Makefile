include ../../.makefile-header

CCFLAGSLOCAL = -I../include

####### Files

SOURCES =   share.c \
            3Dbase.c \
            3Dextra.c \
            2Dbase.c \
            2Dextra.c \
	    1Dbase.c 
OBJECTS =   share.o \
            3Dbase.o \
            3Dextra.o \
            2Dbase.o \
            2Dextra.o \
	    1Dbase.o
TARGETS =   libgridg.a libgridg.so

####### Build rules

all: static
cleanall: makecleanall

static: libgridg.a
shared: libgridg.so

libgridg.a: $(OBJECTS)
	rm -f libgridg.a
	ar cqs libgridg.a $(OBJECTS)
	cp libgridg.a ../.

libgridg.so: $(OBJECTS)
	rm -f libgridg.so
	$(L) -shared $(LFLAGS) -lm -o libgridg.so $(OBJECTS) 

makecleanall: clean
	( cd doc ; make clean )
	( cd eg ; make clean )


# DO NOT DELETE THIS LINE -- make depend depends on it.
