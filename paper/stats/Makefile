PDFLATEX = pdflatex -interaction=nonstopmode

all: main.pdf supplement.pdf

main.pdf: main.tex references.bib
	$(PDFLATEX) main.tex
	bibtex main
	$(PDFLATEX) main.tex
	$(PDFLATEX) main.tex

supplement.pdf: supplement.tex references.bib
	$(PDFLATEX) supplement.tex
	bibtex supplement
	$(PDFLATEX) supplement.tex
	$(PDFLATEX) supplement.tex

.PHONY: all clean simulations
clean:
	rm -f main.aux main.bbl main.blg main.log main.out main.toc \
	      supplement.aux supplement.bbl supplement.blg supplement.log \
	      supplement.out supplement.toc

simulations:
	cd code && python3 run_all.py
