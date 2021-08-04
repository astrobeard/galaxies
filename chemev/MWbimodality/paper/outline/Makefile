
all: outline.pdf 

outline.pdf: outline.tex outline.bib 
	@ pdflatex outline 
	@ bibtex outline 
	@ pdflatex outline 
	@ pdflatex outline 

clean: 
	@ rm -f *.aux 
	@ rm -f *.bbl 
	@ rm -f *.blg 
	@ rm -f *.log 
	@ rm -f *.out 
	@ rm -f *.pdf 
