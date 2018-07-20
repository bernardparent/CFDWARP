include ../../.makefile-header

CCFLAGSLOCAL = -I../include

####### Files

SOURCES =   soap.c printf.c
OBJECTS =   soap.o printf.o
TARGETS =   libsoap.a libsoap.so

####### Build rules

all: static
cleanall: makecleanall

static: libsoap.a
shared: libsoap.so

libsoap.a: $(OBJECTS)
	rm -f libsoap.a
	ar cqs libsoap.a $(OBJECTS)
	cp libsoap.a ../.

libsoap.so: $(OBJECTS)
	rm -f libsoap.so
	$(L) -shared $(LFLAGS) -o libsoap.so $(OBJECTS)

makecleanall: clean
	( cd eg ; make clean )
	( cd doc ; make clean )

# DO NOT DELETE THIS LINE -- make depend depends on it.
