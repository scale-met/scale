TOPDIR      = $(abspath ../..)
include $(TOPDIR)/Mkinclude

all: dox

dox:
	ruby $(TOPDIR)/utils/mkdoc.rb .

clean:
	rm -rf dox

.PHONY: dox clean
