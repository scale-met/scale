###########################################################
#
# Makefile for compiling tex document
#
###########################################################

# target file name
TEXNAME  = scale-gm_users_guide
VERFNAME = 00_version.tex
VERTEMPL = 00_version_template.tex

# commands to compile documents
LATEX    = platex
BIBTEX   = bibtex
DVIPDF   = dvipdfmx -p a4
HTMLBIN  = pandoc --to=html

VERSION  = $(shell cat ../../../scale-rm/src/VERSION)

SRCS = $(wildcard *.tex)

all: pdf

pdf: $(TEXNAME).pdf

html: $(TEXNAME).html

clean:
	rm -f *.pdf *.dvi *.aux *.bbl *.blg *.log *.out *.toc $(VERFNAME)

$(TEXNAME).pdf: $(TEXNAME).dvi
	$(DVIPDF) $(TEXNAME)

$(TEXNAME).dvi: $(VERFNAME) $(SRCS)
	$(LATEX)  $(TEXNAME)
	$(LATEX)  $(TEXNAME)
	$(LATEX)  $(TEXNAME)

$(VERFNAME): $(VERTEMPL)
	$(shell cat $< | sed s/"#VERSION#"/$(VERSION)/g > $@)

$(TEXNAME).html: $(TEXNAME).tex
	$(HTMLBIN) $< -o $@


.SUFFIXES:
.SUFFIXES: .tex .bbl .dvi .pdf .bib .html

.PHONY : pdf clean html
