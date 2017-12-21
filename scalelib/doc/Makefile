################################################################################
#
# Makefile for the reference manual
#
################################################################################

TOPDIR = ../..

include $(TOPDIR)/Mkinclude

all: dox

dox:
	ruby $(TOPDIR)/doc/reference/mkdoc.rb ../src

clean:
	rm -rf dox

.PHONY: dox clean
