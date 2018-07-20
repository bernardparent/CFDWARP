all: src  

src:
	( cd .active ; make src )
	( cd .active ; make arc )

clean:
	( cd .active; make clean )

cleanall:
	( cd _none; make cleanall )
	( cd _explicit; make cleanall )
	( cd _block_ADI; make cleanall )
	( cd _block_DDADI ; make cleanall )
	( cd _block_IMAF ; make cleanall)
	( cd doc ; make clean )

# DO NOT DELETE THIS LINE -- make depend depends on it.
