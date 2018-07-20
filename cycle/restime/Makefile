all: src

src:
	( cd .active ; make src )
	( cd .active ; make arc )

clean:
	( cd .active; make clean )

cleanall:
	( cd _steady-state ; make cleanall )
	( cd _RK4 ; make cleanall )
	( cd _TVDRK ; make cleanall )
	( cd _Trapezoidal ; make cleanall )
	( cd _BDF ; make cleanall )
	( cd _CDF+ ; make cleanall )

# DO NOT DELETE THIS LINE -- make depend depends on it.
