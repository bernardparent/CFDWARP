all: src

src:
	( cd .active ; make src )
	( cd .active ; make arc )

clean:
	( cd .active; make clean )

cleanall:
	( cd doc ; make clean )
	( cd _none ; make cleanall )
	( cd _FVS ; make cleanall )
	( cd _FDS_Cauchy ;  make cleanall )
	( cd _FDS ; make cleanall )
	( cd _FDS+ ; make cleanall )
	( cd _FVSTVD+ ; make cleanall )
	( cd _FVS+ ; make cleanall )
	( cd _MUSCL+ ; make cleanall )
	( cd _MUSCLX ; make cleanall )
	( cd _MUSCL ; make cleanall )

# DO NOT DELETE THIS LINE -- make depend depends on it.
