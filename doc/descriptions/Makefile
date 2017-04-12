###########################################################
#
# Makefile for compiling tex document
#
###########################################################

# commands to compile documents
LATEX    = platex
BIBTEX   = bibtex
DVIPDF   = dvipdfmx -p a4

# target file name
TARGET = scale_rm_description

SRCS = ${TARGET}.tex \
       introduction.tex \
       governing_equation.tex \
       temporal_integration.tex \
       dynamical_process.tex \
       dynamical_process_terrainfollowing.tex \
       dynamical_process_mapfactor.tex \
       dynamical_process_hevi.tex \
       dynamical_process_hivi.tex \
       turbulence.tex \
       turbulence_mynn.tex \
       surface_flux.tex \
       aerosol.tex \
       large_scale_sinking.tex \
       appendix.tex \
       notation.tex \
       variables.tex

all: ${TARGET}.pdf

clean:
	rm -f ${TARGET}.aux ${TARGET}.bbl ${TARGET}.blg ${TARGET}.dvi ${TARGET}.log ${TARGET}.out ${TARGET}.toc

${TARGET}.pdf: ${TARGET}.dvi
	$(DVIPDF) ${TARGET}

${TARGET}.dvi: ${SRCS} reference.bib
	$(LATEX)  ${TARGET}
	$(BIBTEX) ${TARGET}
	$(LATEX)  ${TARGET}
	$(LATEX)  ${TARGET}


