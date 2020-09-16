###########################################################
#
# Makefile for compiling tex document
#
###########################################################

# commands to compile documents
LATEX    = platex
BIBTEX   = bibtex
DVIPDF   = dvipdfmx -p a4

CURL     = curl -f -s -S --create-dirs
IMAGE_URL = https://scale.riken.jp/archives/figure
IMAGE_PATH = $(IMAGE_URL)/$(if $(shell git branch --contains 2>/dev/null|grep -v master|grep -v release|grep -v hotfix),devel,$(VERSION))

VERSION  = $(shell cat ../../scale-rm/src/VERSION)

# target file name
TARGET = scale_rm_description

IMAGE_BASENAME = $(shell grep includegraphics *tex | ruby -n -e '/^\s*\\includegraphics.*\{(.+)\}/ =~ $$_.split(":")[1] && print(/\A([^.]+)(\..+)?\Z/=~File.basename($$1) && $$1+($$2 || ".eps"),"\n") ')
FIG_DIR = ../figure
IMAGES = $(patsubst %,$(FIG_DIR)/%,$(IMAGE_BASENAME))

SRCS = $(wildcard *.tex) $(IMAGES)

all: ${TARGET}.pdf

update_images:
	$(if $(IMAGE_URL),,$(error IMAGE_URL is empty))
	for image in $(IMAGES) ; do \
	  $(CURL) $(IMAGE_PATH)/`basename $$image` -z $$image -o $$image; \
	done

clean:
	rm -f ${TARGET}.aux ${TARGET}.bbl ${TARGET}.blg ${TARGET}.dvi ${TARGET}.log ${TARGET}.out ${TARGET}.toc

${TARGET}.pdf: ${TARGET}.dvi
	$(DVIPDF) ${TARGET}

${TARGET}.dvi: ${SRCS} ../reference.bib
	$(LATEX)  ${TARGET}
	$(BIBTEX) ${TARGET}
	$(LATEX)  ${TARGET}
	$(LATEX)  ${TARGET}

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
