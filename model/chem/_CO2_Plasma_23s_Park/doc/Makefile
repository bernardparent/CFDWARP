####### Build rules

SHELL=/bin/sh

LATEXDIR=".:../../../../lib/latex:"

all: makedvi makebib 
#all: makedvi makebib makeps makepdf

makebib:
	( export BSTINPUTS=$(LATEXDIR) ; exec bibtex report )

makedvi:
	( export TEXINPUTS=$(LATEXDIR) ; exec pdflatex report )

makeps:
	( export TEXINPUTS=$(LATEXDIR) ; exec dvips -t letter report.dvi -o report.ps )

makepdf:
	( exec ps2pdf13 report.ps report.pdf )

clean:
	( rm -f report.dvi )
	( rm -f report.ps )
	( rm -f report.pdf ) 
	( rm -f core )
	( rm -f *~ )
	( rm -f *.aux *.bbl *.blg *.log *.toc *.ilg *.ind) 
	 
     
# DO NOT DELETE THIS LINE -- make depend depends on it.
