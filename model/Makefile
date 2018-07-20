all: src  

src:
	( cd fluid ; make src )
	( cd chem ; make src )
	( cd thermo ; make src )
	( cd metrics ; make src )
	( cd emfield ; make src )
	( cd beam ; make src )
	( cd share ; make src )
	( cd share ; make arc )
	( cd .active ; make src )
	( cd .active ; make arc )

clean:
	( cd .active; make clean )
	( cd fluid ; make clean )
	( cd chem ; make clean )
	( cd share ; make clean )
	( cd thermo ; make clean )
	( cd metrics ; make clean )
	( cd emfield ; make clean )
	( cd beam ; make clean )

cleanall:
	( cd _generic; make cleanall )
	( cd fluid ; make cleanall )
	( cd chem ; make cleanall )
	( cd share ; make cleanall )
	( cd thermo ; make cleanall )
	( cd metrics ; make cleanall )
	( cd emfield ; make cleanall )
	( cd beam ; make cleanall )


# DO NOT DELETE THIS LINE -- make depend depends on it.
