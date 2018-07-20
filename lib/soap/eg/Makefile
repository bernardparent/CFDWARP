include ../../../.makefile-header

CCFLAGSLOCAL = -I../../include
LFLAGSLOCAL = -L../../ -lsoap -lexm -lm -lpthread

####### Files

SOURCES =   bare.c variables.c functions.c functions2.c actions.c array.c calc.c soap.c
OBJECTS =   bare.o variables.o functions.o functions2.o actions.o array.o calc.o soap.o
TARGETS =   bare variables functions functions2 actions array calc soap

####### Build rules


all: $(TARGETS)

bare: bare.o
	$(L) bare.o $(LFLAGS) $(LFLAGSLOCAL) -o bare

variables: variables.o
	$(L) variables.o $(LFLAGS) $(LFLAGSLOCAL) -o variables

functions: functions.o
	$(L) functions.o $(LFLAGS) $(LFLAGSLOCAL) -o functions

functions2: functions2.o
	$(L) functions2.o $(LFLAGS) $(LFLAGSLOCAL) -o functions2

actions: actions.o
	$(L) actions.o $(LFLAGS) $(LFLAGSLOCAL) -o actions

calc: calc.o
	$(L) calc.o $(LFLAGS) $(LFLAGSLOCAL) -o calc

array: array.o
	$(L) array.o $(LFLAGS) $(LFLAGSLOCAL) -o array

soap: soap.o
	$(L) soap.o $(LFLAGS) $(LFLAGSLOCAL) -o soap

# DO NOT DELETE THIS LINE -- make depend depends on it.
