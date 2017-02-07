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
	$(MAKE) -C $(TOPDIR)/scale-rm/doc dox

doxygen:
	doxygen scale-rm.cfg

clean:
	rm -rf html latex

allclean: clean
	$(MAKE) -C $(TOPDIR)/scalelib/doc clean

.PHONY: doxygen clean
