PDFLATEX = pdflatex -interaction=nonstopmode

main.pdf: main.tex references.bib
	$(PDFLATEX) main.tex
	bibtex main
	$(PDFLATEX) main.tex
	$(PDFLATEX) main.tex

.PHONY: clean simulations
clean:
	rm -f main.aux main.bbl main.blg main.log main.out main.toc

simulations:
	cd code && python3 run_all.py
