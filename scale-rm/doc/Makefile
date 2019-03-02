################################################################################
#
# Makefile for the reference manual
#
################################################################################

TOPDIR = ../..

all: dox

dox:
	ruby $(TOPDIR)/doc/reference/mkdoc.rb ../src

clean:
	rm -rf dox

.PHONY: dox clean
