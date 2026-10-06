# Run make in this directory. KUDPC users should run the payload on a compute node.
XELATEX ?= xelatex
BIBTEX ?= bibtex
DOCUMENT := beach_models
BUILD := build
SOURCES := $(DOCUMENT).tex references.bib $(wildcard sections/*.tex)

.DELETE_ON_ERROR:
.PHONY: all clean
all: $(BUILD)/$(DOCUMENT).pdf

$(BUILD)/$(DOCUMENT).pdf: $(SOURCES)
	mkdir -p $(BUILD)
	$(XELATEX) -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=$(BUILD) $(DOCUMENT).tex
	cd $(BUILD) && BIBINPUTS=..: $(BIBTEX) $(DOCUMENT)
	$(XELATEX) -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=$(BUILD) $(DOCUMENT).tex
	$(XELATEX) -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=$(BUILD) $(DOCUMENT).tex

clean:
	rm -f $(BUILD)/$(DOCUMENT).aux $(BUILD)/$(DOCUMENT).bbl $(BUILD)/$(DOCUMENT).blg
	rm -f $(BUILD)/$(DOCUMENT).log $(BUILD)/$(DOCUMENT).out $(BUILD)/$(DOCUMENT).toc
	rm -f $(BUILD)/$(DOCUMENT).pdf
