NUM ?= 3 # For LaTeX to correctly cross-reference.  To save time, use 1.
VERBOSE ?= 0
ifeq ($(VERBOSE),0)
	CMDLOG_DIR_REDIRECT := > /dev/null
	PSTAKE := ./pstake.py -q
else
	CMDLOG_DIR_REDIRECT := >&1
	PSTAKE := ./pstake.py
endif

# tool binaries
LATEX := latex
BIBTEX := bibtex
DVIPS := dvips
PS2PDF := ps2pdf

# directories
LOG_DIR := log
EPS_DIR := thingol_eps
SCHEMATIC_DIR := schematic

TRUNK := thingol

MAKEFILES := Makefile

ALL_TEX = $(wildcard $(SCHEMATIC_DIR)/*.tex)
ALL_EPS = $(patsubst $(SCHEMATIC_DIR)/%.tex,$(EPS_DIR)/%.eps,$(ALL_TEX))

.PHONY: default
default: $(TRUNK)_mesh.pdf

$(EPS_DIR)/%.eps: $(SCHEMATIC_DIR)/%.tex $(MAKEFILES) pstake.py
	mkdir -p $(EPS_DIR)
	$(PSTAKE) $< $@

$(TRUNK)_mesh.dvi: $(TRUNK)_mesh.tex ../$(TRUNK)_main.bib $(ALL_EPS) $(MAKEFILES)
	mkdir -p $(LOG_DIR)
	@echo "Having EPS files: $(ALL_EPS)"
	num=1 ; while [ $$num -le $(NUM) ] ; do \
		$(LATEX) $< 2>&1 | tee $(LOG_DIR)/$@.latex.$$num.log $(CMDLOG_DIR_REDIRECT) ; \
		$(BIBTEX) ../$(basename $<) 2>&1 | tee $(LOG_DIR)/$@.bibtex.$$num.log $(CMDLOG_DIR_REDIRECT) ; \
	(( num = num + 1 )) ; done

$(TRUNK)_mesh.ps: $(TRUNK)_mesh.dvi $(MAKEFILES)
	mkdir -p $(LOG_DIR)
	$(DVIPS) $< 2>&1 | tee $(LOG_DIR)/$@.log $(CMDLOG_DIR_REDIRECT)

$(TRUNK)_mesh.pdf: $(TRUNK)_mesh.ps $(MAKEFILES)
	mkdir -p $(LOG_DIR)
	$(PS2PDF) $< 2>&1 | tee $(LOG_DIR)/$@.log $(CMDLOG_DIR_REDIRECT)

.PHONY: clean_tex
clean_tex:
	rm -f \
		$(TRUNK)_*.aux \
		$(TRUNK)_*.bbl \
		$(TRUNK)_*.blg \
		$(TRUNK)_*.dvi \
		$(TRUNK)_*.log \
		$(TRUNK)_*.out \
		$(TRUNK)_*.pdf \
		$(TRUNK)_*.ps

.PHONY: clean_eps
clean_eps:
	rm -rf $(LOG_DIR) $(EPS_DIR)

.PHONY: clean
clean: clean_tex clean_eps
