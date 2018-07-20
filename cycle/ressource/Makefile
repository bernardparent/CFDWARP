all: src

src:
	( cd .active ; make src )
	( cd .active ; make arc )

clean:
	( cd .active; make clean )

cleanall:
	( cd _forward ; make cleanall )
	( cd _SVS+ ; make cleanall )
	( cd _SDS+ ; make cleanall )
	( cd doc ; make clean )

# DO NOT DELETE THIS LINE -- make depend depends on it.
