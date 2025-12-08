FIGURES=$(patsubst %.svg,%.pdf,$(wildcard figures/*.svg))
SVG_RESOLUTION=-density 300

all: Slides

figures/%.pdf: figures/%.svg
	convert $(SVG_RESOLUTION) $< $@


Slides: Slides.tex $(FIGURES)
	pdflatex Slides
	pdflatex Slides

.phony: clean
clean:
	rm -f ${FIGURES}