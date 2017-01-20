TARGET = scale_rm_description

all: ${TARGET}.pdf

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

${TARGET}.pdf: ${SRCS} ${TARGET}.aux ${TARGET}.bib
	latex $<
	dvipdfmx ${TARGET}

${TARGET}.aux: ${SRCS}
	latex ${TARGET}

${TARGET}.bib: ${SRCS} ${TARGET}.aux
	bibtex ${TARGET}


clean:
	rm -f ${TARGET}.aux ${TARGET}.bbl ${TARGET}.blg ${TARGET}.dvi ${TARGET}.log ${TARGET}.out ${TARGET}.toc
