################################################################################
#
# Makefile for the reference manual
#
################################################################################

TOPDIR      = $(abspath ../..)

include $(TOPDIR)/Mkinclude

all: dox doxygen tex

dox:
	$(MAKE) -C $(TOPDIR)/scalelib/doc dox
	$(MAKE) -C $(TOPDIR)/scale-rm/doc dox

doxygen:
	doxygen scale-rm.cfg

tex: latex/refman.pdf

latex/refman.pdf: latex/refman.tex
	$(MAKE) -C latex

clean:
	rm -rf html latex

allclean: clean
	$(MAKE) -C $(TOPDIR)/scalelib/doc clean

.PHONY: doxygen clean
