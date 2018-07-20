####### Build rules

default: 
	( echo 'type make cleanout to clean up this directory from outfiles and postfiles' )

cleanout: makecleanout

makecleanout: 
	( rm -f tecplot.phy `find . -type f -name 'tecplot.phy' -print` ) 
	( rm -f gmon.out `find . -type f -name 'gmon.out' -print` ) 
	( rm -f batch.log `find . -type f -name 'batch.log' -print` ) 
	( rm -f tp* `find . -type f -name 'tp*' -print` ) 
	( rm -f out* `find . -type f -name 'out.*' -print` ) 
	( rm -f post* `find . -type f -name 'post.*' -print` ) 
	( rm -f post* `find . -type f -name 'gprof.*' -print` ) 
	( rm -f scr* `find . -type f -name 'scr.*' -print` ) 
	( rm -f *ps `find . -type f -name '*ps' -print` ) 
	( rm -f *.avi `find . -type f -name '*.avi' -print` ) 
	( rm -f *.mp4 `find . -type f -name '*.mp4' -print` ) 
	( rm -f *.mpg `find . -type f -name '*.mpg' -print` ) 

clean:
	( rm -f tecplot.phy `find . -type f -name 'tecplot.phy' -print` ) 
	 
     
# DO NOT DELETE THIS LINE -- make depend depends on it.
