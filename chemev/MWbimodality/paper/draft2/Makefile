
all: ms.pdf 

ms.pdf: *.tex ms.bib 
	pdflatex ms 
	bibtex ms 
	pdflatex ms 
	pdflatex ms 

.PHONY: clean 
clean: 
	@ rm -f *.aux 
	@ rm -f *.bcf 
	@ rm -f *.blg 
	@ rm -f *.log 
	@ rm -f *.out 
	@ rm -f *.run.xml 
	@ rm -f ms.pdf 

