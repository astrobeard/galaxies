
all: outline2.pdf 

outline2.pdf: outline2.tex outline2.bib 
	@ pdflatex outline2 
	@ bibtex outline2 
	@ pdflatex outline2 
	@ pdflatex outline2 

clean: 
	@ rm -f *.aux 
	@ rm -f *.bbl 
	@ rm -f *.blg 
	@ rm -f *.log 
	@ rm -f *.out 
	@ rm -f *.pdf 
