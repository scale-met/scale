TOPDIR = ../..
include $(TOPDIR)/Mkinclude

all: dox doxygen


dox:
	$(MAKE) -C ../../scalelib/doc dox
	ruby $(TOPDIR)/utils/mkdoc.rb .

doxygen:
	doxygen scaleles.cfg

clean:
	rm -rf dox html latex

allclean: clean
	$(MAKE) -C ../../scalelib/doc clean

.PHONY: dox doxygen clean
