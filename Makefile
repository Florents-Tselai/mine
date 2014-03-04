VIRTUAL_ENV_DIR=mine_env
GPROF2DOT_PATH=mine_env/lib/python2.7/site-packages/gprof2dot
PROFILING_RESULTS_DIR=doc/profiling
EXPERIMENTS_DIR=doc/experiments
EXAMPLES_DIR=doc/examples
create_virtualenv:
	virtualenv $(VIRTUAL_ENV_DIR) && \
	. $(VIRTUAL_ENV_DIR)/bin/activate && \
	pip install matplotlib numpy gprof2dot pandas

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
	mkdir -p $(PROFILING_RESULTS_DIR)
	python -m cProfile --sort cumulative -o mine.pstats mine/test_mine.py
	python $(GPROF2DOT_PATH)/gprof2dot.py -f pstats mine.pstats | dot -Tpdf -o $(PROFILING_RESULTS_DIR)/profile.pdf

test:
	python mine/test_mine.py

run_experiments:
	mkdir -p $(EXPERIMENTS_DIR)
	python mine/experiments.py

plot_examples:
	mkdir -p $(EXAMPLES_DIR)
	python mine/examples.py
