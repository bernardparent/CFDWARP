all: src

src:
	( cd .active ; make src )
	( cd .active ; make arc )

clean:
	( cd .active ; make clean )

cleanall:
	( cd _none ; make cleanall )
	( cd _ebeam_fixed ; make cleanall )
	( cd _ebeam_algebraic ; make cleanall )
	( cd doc ; make clean )

# DO NOT DELETE THIS LINE -- make depend depends on it.
