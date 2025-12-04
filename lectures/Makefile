FIGURES=$(patsubst %.svg,%.pdf,$(wildcard figures/*.svg))
SVG_RESOLUTION=-density 300

all: Notes

figures/%.pdf: figures/%.svg
	convert $(SVG_RESOLUTION) $< $@


Notes: Notes.tex $(FIGURES)
	pdflatex Notes

.phony: clean
clean:
	rm -f ${FIGURES}