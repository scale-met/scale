###########################################################
#
# Makefile for compiling tex document
#
###########################################################

# target file name
TEXNAME  = scale_users_guide
VERFNAME = 00_version.tex

# commands to compile documents
LATEX    = platex
BIBTEX   = bibtex
DVIPDF   = dvipdfmx -p a4
HTMLBIN  = pandoc --to=html

VERSION  = $(shell cat ../../../src/VERSION)

all: pdf

pdf: $(TEXNAME).pdf

clean:
	rm -f *.pdf *.dvi *.aux *.bbl *.blg *.log *.out *.toc

.dvi.pdf:
	$(DVIPDF) $(TEXNAME)

.tex.dvi:
	$(shell echo "\newcommand{\version}{$(VERSION)}" > $(VERFNAME))
	$(LATEX) $(TEXNAME)
	$(BIBTEX) $(TEXNAME)
	$(LATEX) $(TEXNAME)
	$(LATEX) $(TEXNAME)

html: $(TEXNAME).html
	$(HTMLBIN) $(TEXNAME).tex -o ${TEXNAME}.html

.SUFFIXES:
.SUFFIXES: .tex .bbl .dvi .pdf .bib .html

.PHONY : pdf clean html
