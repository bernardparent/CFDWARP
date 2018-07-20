all: src

src:
	( cd ts ; make src )
	( cd tsemf ; make src )
	( cd res ; make src )
	( cd resconv ; make src )
	( cd restime ; make src )
	( cd ressource ; make src )
	( cd share ; make src )
	( cd share ; make arc )
	( cd .active ; make src )
	( cd .active ; make arc )

clean:
	( cd .active; make clean )
	( cd ts ; make clean )
	( cd tsemf ; make clean )
	( cd res ; make clean )
	( cd restime ; make clean )
	( cd ressource ; make clean )
	( cd resconv ; make clean )
	( cd share ; make clean )

cleanall:
	( cd _multizone_marching; make cleanall )
	( cd _standard; make cleanall )
	( cd _predictor_corrector; make cleanall )
	( cd _multizone; make cleanall )
	( cd ts ; make cleanall )
	( cd tsemf ; make cleanall )
	( cd res ; make cleanall )
	( cd resconv ; make cleanall )
	( cd restime ; make cleanall )
	( cd ressource ; make cleanall )
	( cd share ; make cleanall )
	( cd doc ; make clean )


# DO NOT DELETE THIS LINE -- make depend depends on it.
