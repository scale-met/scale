################################################################################
#
# Makefile for main program
#
################################################################################

TOPDIR      = $(abspath ../..)

include $(TOPDIR)/Mkinclude

all: dox doxygen

dox:
	$(MAKE) -C $(TOPDIR)/scalelib/doc dox
	ruby $(TOPDIR)/utils/mkdoc.rb .

doxygen:
	doxygen scale-rm.cfg

clean:
	rm -rf dox html latex

allclean: clean
	$(MAKE) -C $(TOPDIR)/scalelib/doc clean

.PHONY: dox doxygen clean
