include ../../.makefile-header

CCFLAGSLOCAL = -I../include

####### Files

SOURCES =   exm.c
OBJECTS =   exm.o
TARGETS =   libexm.a libexm.so

####### Build rules

all: static

static: libexm.a
shared: libexm.so

libexm.a: $(OBJECTS)
	rm -f libexm.a
	ar cqs libexm.a $(OBJECTS)
	cp libexm.a ../.

libexm.so: $(OBJECTS)
	rm -f libexm.so
	$(L) -shared $(LFLAGS) -o libexm.so $(OBJECTS) 

# DO NOT DELETE THIS LINE -- make depend depends on it.
