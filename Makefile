doc:
	cd doc;\
	latex mic.tex;\
	bibtex mic;\
	pdflatex mic.tex;\
	pdflatex mic.tex
	
