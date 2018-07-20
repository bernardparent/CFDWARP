####### Build rules

SHELL=/bin/sh

LATEXDIR=".:../../../../lib/latex:"

all: makedvi makebib makeps

makebib:
	( export BIBINPUTS=$(LATEXDIR) ; exec bibtex report )

makedvi:
	( export TEXINPUTS=$(LATEXDIR) ; exec latex report )

makeps:
	( export TEXINPUTS=$(LATEXDIR) ; exec dvips -t letter report.dvi -o report.ps )

clean:
	( rm -f report.dvi )
	( rm -f report.ps )
	( rm -f core )
	( rm -f *.aux *.bbl *.blg *.log *.toc *.ilg *.ind) 
	 
     
# DO NOT DELETE THIS LINE -- make depend depends on it.
