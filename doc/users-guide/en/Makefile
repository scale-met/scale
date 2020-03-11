###########################################################
#
# Makefile for compiling tex document
#
###########################################################

# target file name
TEXNAME  = scale_users_guide_En
VERFNAME = 00_version.tex
VERTEMPL = 00_version_template.tex

# commands to compile documents
LATEX    = pdflatex
BIBTEX   = bibtex
HTMLBIN  = pandoc --to=html
CURL     = curl -f -s -S --create-dirs
IMAGE_URL = http://scale.aics.riken.jp/doc/figure
IMAGE_PATH = $(IMAGE_URL)/$(if $(shell git branch --contains 2>/dev/null|grep -v master|grep -v release|grep -v hotfix),devel,$(VERSION))

VERSION  = $(shell cat ../../../scale-rm/src/VERSION)

IMAGE_BASENAME = $(shell grep includegraphics *tex | ruby -n -e '/^\s*\\includegraphics.*\{(.+)\}/ =~ $$_.split(":")[1] && print(/\A([^.]+)(\..+)?\Z/=~File.basename($$1) && $$1+($$2 || ".eps"),"\n") ') title_wallpaper.pdf scale_logo_final_ULWB.pdf
FIG_DIR = ../../figure
IMAGES = $(patsubst %,$(FIG_DIR)/%,$(IMAGE_BASENAME))

SRCS = $(wildcard *.tex) $(IMAGES)


all: pdf

pdf: $(TEXNAME).pdf

html: $(TEXNAME).html

update_images:
	$(if $(IMAGE_URL),,$(error IMAGE_URL is empty))
	for image in $(IMAGES) ; do \
	  $(CURL) $(IMAGE_PATH)/`basename $$image` -z $$image -o $$image; \
	done

clean:
	rm -f *.dvi *.aux *.bbl *.blg *.log *.out *.toc $(VERFNAME)

allclean: clean
	rm -rf *.pdf $(FIG_DIR)

$(TEXNAME).pdf: $(VERFNAME) $(SRCS) ../../reference.bib
	$(LATEX)  $(TEXNAME)
	$(BIBTEX) $(TEXNAME)
	$(LATEX)  $(TEXNAME)
	$(LATEX)  $(TEXNAME)

$(VERFNAME): $(VERTEMPL)
	$(shell cat $< | sed s/"#VERSION#"/$(VERSION)/g > $@)

$(TEXNAME).html: $(TEXNAME).tex
	$(HTMLBIN) $< -o $@


.SUFFIXES:
.SUFFIXES: .tex .bbl .dvi .pdf .bib .html .png .eps

$(FIG_DIR)/%.png:
	$(if $(IMAGE_URL),,$(error IMAGE_URL is empty))
	$(CURL) $(IMAGE_PATH)/`basename $@` -o $@
$(FIG_DIR)/%.pdf:
	$(if $(IMAGE_URL),,$(error IMAGE_URL is empty))
	$(CURL) $(IMAGE_PATH)/`basename $@` -o $@
$(FIG_DIR)/%.eps:
	$(if $(IMAGE_URL),,$(error IMAGE_URL is empty))
	$(CURL) $(IMAGE_PATH)/`basename $@` -o $@

.PHONY : pdf clean html
