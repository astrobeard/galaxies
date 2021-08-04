
all: figures.pdf 

%.pdf: %.tex 
	@ pdflatex $< 
	@ pdflatex $< 

.PHONY: clean 
clean: 
	@ rm -f *.aux 
	@ rm -f *.log 
	@ rm -f *.pdf 

