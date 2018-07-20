####### Build rules

all: makelib makesrc maketar maketools makedoc 

lib: makelib
config: makeconfig
configx: makeconfigx
configp: makeconfigp
configpx: makeconfigpx
src: makelib makesrc
tools: makelib maketools
tar: makelib makesrc maketar

makelib:
	( cd lib ; make all )

makelibext:
	( cd lib ; make ext )

makeconfig:
	( cd config ; make clean ; make config )
	( cd config ; exec ./config )
	( cd config ; chmod u+x applyconfig.sh )
	( cd config ; exec ./applyconfig.sh )

makeconfigx:
	( cd config ; make clean ; make configx )
	( cd config ;  exec ./config )
	( cd config ; chmod u+x applyconfig.sh )
	( cd config ; exec ./applyconfig.sh )

makeconfigp:
	( cd config ; make clean ; make configp )
	( cd config ;  exec ./config )
	( cd config ; chmod u+x applyconfig.sh )
	( cd config ; exec ./applyconfig.sh )

makeconfigpx:
	( cd config ; make clean ; make configpx )
	( cd config ;  exec ./config )
	( cd config ; chmod u+x applyconfig.sh )
	( cd config ; exec ./applyconfig.sh )

makesrc:	
	( cd src ; make )
	@echo
	@echo Finished compiling CFDWARP.
	@echo 
	@echo The execuble is named warp and is located in the src directory.
	@echo Type warp -h for a list of command line options.
	@echo 

maketar:	
	( cd tar ; make )

maketools:	
	( cd tools ; make ) 

cleanbin:	
	( cd bin ; make cleanout ) 

clean:
	( cd config ; make clean ) 
	( cd src ; make clean ) 
	( cd model ; make clean ) 
	( cd cycle ; make clean )
	( cd tools ; make clean ) 
	( cd lib ; make clean ) 
	( cd tar ; make clean ) 
	( rm -f core `find . -type f -name 'core' -print` ) 
	@echo
	@echo Finished cleaning the directories.
	@echo 
	@echo Type make src to compile CFDWARP.
	@echo 


cleanall:
	( cd config ; make cleanall ) 
	( cd src ; make cleanall ) 
	( cd tar ; make cleanall ) 
	( cd model ; make cleanall ) 
	( cd cycle ; make cleanall ) 
	( cd tools ; make cleanall ) 
	( cd lib ; make cleanall )
	( rm -f core `find . -type f -name 'core' -print` )
	( rm -f gmon.out `find . -type f -name 'gmon.out' -print` )
	( rm -f  `find . -type f -name '*~' -print` )
	( cp .makefile-header-default .makefile-header )


# DO NOT DELETE THIS LINE -- make depend depends on it.
