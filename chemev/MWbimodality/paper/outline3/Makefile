
all: outline3.pdf 

outline3.pdf: outline3.tex outline3.bib 
	@ pdflatex outline3 
	@ bibtex outline3 
	@ pdflatex outline3 
	@ pdflatex outline3 

clean: 
	@ rm -f *-blx.bib 
	@ rm -f *.aux 
	@ rm -f *.bbl 
	@ rm -f *.blg 
	@ rm -f *.log 
	@ rm -f *.out 
	@ rm -f *.pdf 
