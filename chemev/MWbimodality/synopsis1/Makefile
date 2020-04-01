
COMPILER = pdflatex 

all: synopsis.pdf 
.PHONY: clean 

%.pdf: %.tex 
	@ $(COMPILER) $< 
	@ $(COMPILER) $< 

clean: 
	@ rm -f *.aux 
	@ rm -f *.log 
	@ rm -f *.out 
	@ rm -f *.pdf 

