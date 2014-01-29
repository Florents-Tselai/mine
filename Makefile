doc:
	cd pub;\
	latex mic.tex;\
	bibtex mic;\
	pdflatex mic.tex;\
	pdflatex mic.tex
	
