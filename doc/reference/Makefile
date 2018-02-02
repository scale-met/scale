################################################################################
#
# Makefile for the reference manual
#
################################################################################

TOPDIR = ../..

include $(TOPDIR)/Mkinclude

all: doxygen tex

dox:
	$(MAKE) -C $(TOPDIR)/dc_utils/doc dox
	$(MAKE) -C $(TOPDIR)/scalelib/doc dox
	$(MAKE) -C $(TOPDIR)/scale-rm/doc dox

doxygen: dox
	doxygen scale-rm.cfg

tex: latex/refman.pdf

latex/refman.pdf: latex/refman.tex
	$(MAKE) -C latex

clean:
	rm -rf html latex

allclean: clean
	$(MAKE) -C $(TOPDIR)/dc_utils/doc clean
	$(MAKE) -C $(TOPDIR)/scalelib/doc clean
	$(MAKE) -C $(TOPDIR)/scale-rm/doc clean

.PHONY: doxygen clean
