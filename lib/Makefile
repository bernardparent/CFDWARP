all: makeexm makesoap makegridg

exm: makeexm
soap: makesoap
gridg: makeexm makesoap makegridg

makeexm:
	( cd exm ; make )

makesoap:
	( cd soap ; make )

makegridg:
	( cd gridg ; make )

clean:
	( cd gridg; make clean )
	( cd exm; make clean )
	( cd soap; make clean )
	( rm -f *.a )

cleanall:
	( cd gridg; make cleanall )
	( cd exm; make cleanall )
	( cd soap; make cleanall )
	( cd latex/eg; make clean )
	( rm -f *.a )

# DO NOT DELETE THIS LINE -- make depend depends on it.
