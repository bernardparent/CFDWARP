SHELL=/bin/sh

LATEXDIR=".:../../lib/latex:"                                                                                    

all: makedvi makeps makebib makepdf

makebib:
	( export BSTINPUTS=$(LATEXDIR) ; export BSTINPUTS=$(LATEXDIR) ; exec bibtex report )

makedvi:
	( export TEXINPUTS=$(LATEXDIR) ; exec latex report )

makeps:
	( export TEXINPUTS=$(LATEXDIR) ; exec dvips -t letter report.dvi -o report.ps )

makepdf:
	( exec ps2pdf13 report.ps report.pdf )

clean:
	( rm -f core missfont.log *.bak *~ )
	( rm -f report.dvi report.pdf report.ps report.aux report.blg report.log report.toc report.ttt report.fff report.lot report.lof )


