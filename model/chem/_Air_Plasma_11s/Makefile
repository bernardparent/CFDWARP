include ../../../.makefile-header

LFLAGSLOCAL =  
CCFLAGSLOCAL = $(CFLAGSCFD) -I../../../lib/include -I../../../.

####### Files

HEADERS =   chem.h dunnkang1973.h park1993.h boyd2007.h lenard1964.h farbar2013.h parentdunn2021.h  parentpark2021.h kim2021.h parent2023.h parent2023b.h thoguluva2023.h rodriguez2024.h rodriguezpark2024.h rodriguez2026.h
SOURCES =   chem.c dunnkang1973.c park1993.c boyd2007.c lenard1964.c farbar2013.c parentdunn2021.c  parentpark2021.c kim2021.c parent2023.c parent2023b.c thoguluva2023.c rodriguez2024.c rodriguezpark2024.c rodriguez2026.c
OBJECTS =   chem.o dunnkang1973.o park1993.o boyd2007.o lenard1964.o farbar2013.o parentdunn2021.o  parentpark2021.o kim2021.o parent2023.o parent2023b.o thoguluva2023.o rodriguez2024.o rodriguezpark2024.o rodriguez2026.o
TARGETS  = ../_chem.a


####### Build rules

all: src 

src: $(OBJECTS) 

arc: 
	( rm -f $(TARGETS) ; ar cqs $(TARGETS) *.o )                        

cleanall:
	( cd doc ; make clean )
	( make clean )

	

# DO NOT DELETE THIS LINE -- make depend depends on it.
