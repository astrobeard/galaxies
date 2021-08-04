
all: draft1.pdf 

draft1.pdf: draft1.tex draft1.bib 
	pdflatex draft1 
	bibtex draft1 
	pdflatex draft1 
	pdflatex draft1 

.PHONY: clean 
clean: 
	@ rm -f *.aux 
	@ rm -f *.bbl 
	@ rm -f *.bcf 
	@ rm -f *.blg 
	@ rm -f *.log 
	@ rm -f *.out 
	@ rm -f *.run.xml 
	@ rm -f *.pdf 

