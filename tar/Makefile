####### Build rules

SHELL=/bin/sh


all: gridg soap warp tgz

warp:
	( cp ../src/warp warp_package/warp ; cp ../LICENSE warp_package/LICENSE.txt ; cp ../USAGE warp_package/USAGE.txt)

gridg:
	( cd ../lib/gridg/doc/ ; make ; make ; make )
	( mkdir warp_package ; cp ../lib/gridg/doc/report.pdf warp_package/GRIDG_manual.pdf )

soap:
	( cd ../lib/soap/doc/ ; make ; make ; make )
	( mkdir warp_package ; cp ../lib/soap/doc/report.pdf warp_package/SOAP_manual.pdf )

tgz:
	( tar -cpPz warp_package > warp_package.tgz )

clean:
	( rm -f warp_package/GRIDG_manual.pdf ) 
	( rm -f warp_package/SOAP_manual.pdf ) 
	( rm -f warp_package/LICENSE.txt ) 
	( rm -f warp_package/USAGE.txt ) 
	( rm -f warp_package/warp ) 
	( rm -f warp_package.tgz ) 

cleanall: clean 
	 
     
# DO NOT DELETE THIS LINE -- make depend depends on it.
