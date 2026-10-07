all: src

src:
	( cd .active ; make src )
	( cd .active ; make arc )


clean:
	( cd .active ; make clean )

cleanall:
	( cd _Navier-Stokes_perfect ; make cleanall )
	( cd _Navier-Stokes ; make cleanall )
	( cd _Navier-Stokes_plasma ; make cleanall )
	( cd _Favre-Reynolds ; make cleanall )
	( cd _Favre-Reynolds_Tv ; make cleanall )
	( cd _Favre-Reynolds_plasma ; make cleanall )
	( cd _drift-diffusion ; make cleanall )
	( cd doc ; make clean )

# DO NOT DELETE THIS LINE -- make depend depends on it.
