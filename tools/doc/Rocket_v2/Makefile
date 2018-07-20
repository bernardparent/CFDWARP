####### Build rules

SHELL=/bin/sh

LATEXDIR=".:../../../lib/latex:"

all: makedvi makebib makeps makepdf

makebib:
	( export BSTINPUTS=$(LATEXDIR) ; exec bibtex Rocketrep )

makedvi:
	( export TEXINPUTS=$(LATEXDIR) ; exec latex Rocketrep )

makeps:
	( export TEXINPUTS=$(LATEXDIR) ; exec dvips -t letter Rocketrep.dvi -o Rocketrep.ps )

makepdf:
	( ps2pdf13 Rocketrep.ps Rocketrep.pdf )

clean:
	( rm -f Rocketrep.dvi )
	( rm -f Rocketrep.ps )
	( rm -f core )
	( rm -f *.aux *.bbl *.blg *~ *.log *.toc *.ilg *.ind) 
	 
     
# DO NOT DELETE THIS LINE -- make depend depends on it.
