pdf:
	cd doc && \
	pdflatex mic.tex && \
	bibtex mic && \
	pdflatex mic.tex && \
	pdflatex mic.tex

clean:
	cd doc && \
	rm -f mic.aux && \
	rm -f mic.log && \
	rm -f mic.pdf && \
	rm -f mic.bbl && \
	rm -f mic.out && \
	rm -f mic.blg
