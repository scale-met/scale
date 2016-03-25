TOPDIR = ../..
include $(TOPDIR)/Mkinclude

all: list doxygen


list:
	ruby $(TOPDIR)/utils/mkdoc.rb .

doxygen:
	doxygen scaleles.cfg

.PHONY: list doxygen
