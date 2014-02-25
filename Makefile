create_virtualenv:
	virtualenv mine_env && \
	. mine_env/bin/activate
	pip install matplotlib numpy gprof2dot

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

profile:
	python -m cProfile --sort cumulative -o mine.pstats mine/test_mine.py
	python mine_env/lib/python2.7/site-packages/gprof2dot/gprof2dot.py -f pstats mine.pstats | dot -Tpng -o profile.png

