PDFS = design.pdf modulestructure.pdf using-git.pdf

all:	$(PDFS)

install: all
	mkdir -p $(INSTALLPATH)
	cp $(PDFS) $(INSTALLPATH)

%.pdf: %.tex Makefile
	pdflatex $*
	-bibtex $*
	pdflatex $*
	pdflatex $*

%.pdf: %.lyx Makefile
	lyx --export pdf2 $<

%.tex: %.lyx Makefile
	lyx --export pdflatex $<

texclean:	
	rm -f *.aux *.bbl *.blg *.log

clean:	texclean
	rm -f $(PDFS)

.PHONY: all texclean clean
