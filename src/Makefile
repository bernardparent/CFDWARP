include ../.makefile-header

GIT := $(shell command -v git 2> /dev/null)

LFLAGSLOCAL = -L../lib -lgridg -lsoap -lexm -lpthread  -lm

#LFLAGSLOCAL = -L../lib -lgridg -lsoap -lexm -lpthread  -lm -lgomp
#LFLAGSLOCAL = -lm -L../lib -lgridg -lsoap -lexm -lpthread -lnjamd
#LFLAGSLOCAL = -lm -L../lib -lgridg -lsoap -lexm -lpthread -lmemwatch
#LFLAGSLOCAL = -lm -L../lib -lgridg -lsoap -lexm -lpthread -lefence
#LFLAGSLOCAL = -lm -L../lib -lgridg -lsoap -lexm -lpthread -lccmalloc -ldl
#LFLAGSLOCAL = -lm -L../lib -lgridg -lsoap -lexm -lpthread -lmpatrol -lbfd -liberty
#LFLAGSLOCAL = -lm -L../lib -lgridg -lsoap -lexm -lpthread -lmss


CCFLAGSLOCAL = $(CFLAGSCFD) -I../. -I../lib/include

####### Files

HEADERS =   common.h bdry.h control.h post.h data.h init.h
SOURCES =   common.c bdry.c control.c post.c data.c init.c 
OBJECTS =   common.o bdry.o control.o post.o data.o init.o 
TARGETS  =	warp test

####### Build rules

all: rmtargets findversion model cycle warp test

rmtargets:
	( rm -f $(TARGETS) )

findversion:
ifdef GIT
	( echo '#define VERSION "'`git describe --tags`'"' > version.hh )
endif

model:
	( cd ../model ; make src )

cycle:
	( cd ../cycle ; make src )

warp: $(OBJECTS) main.o
	( cd arc ; rm -f *.o ; ar -x ../../model/_model.a ; ar -x ../../cycle/_cycle.a )
	$(L) $(OBJECTS) main.o arc/*.o $(LFLAGS) $(LFLAGSLOCAL) -o warp

test: $(OBJECTS) test.o
	( cd arc ; rm -f *.o ; ar -x ../../model/_model.a ; ar -x ../../cycle/_cycle.a )
	$(L) $(OBJECTS) test.o arc/*.o $(LFLAGS) $(LFLAGSLOCAL) -o test

cleanall:
	( rm -f arc/* )

# DO NOT DELETE THIS LINE -- make depend depends on it.
