####### Build rules

SHELL=/bin/sh

LATEXDIR=".:../../latex:"

all: makedvi makeps makepdf

makedvi:
	( export TEXINPUTS=$(LATEXDIR) ; exec latex report )

makeps:
	( export TEXINPUTS=$(LATEXDIR) ; exec dvips -t letter report.dvi -o report.ps )

makepdf:
	( ps2pdf13 report.ps report.pdf )

clean:
	( rm -f report.dvi ) 
	( rm -f report.ps ) 
	( rm -f report.pdf ) 
	( rm -f core ) 
	( rm -f *.aux *.bbl *~ *.blg *.log *.toc *.bak ) 
	 
     
# DO NOT DELETE THIS LINE -- make depend depends on it.
