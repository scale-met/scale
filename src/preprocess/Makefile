SUBDIRS=boundary grid initial refstate

all:
	@for i in $(SUBDIRS) ; do \
	  cd $$i ; \
	  echo "Making $$i/" ; \
	  $(MAKE) ; \
	  cd .. ; \
	done

clean:
	@for i in $(SUBDIRS) ; do \
	  cd $$i ; \
	  echo "Cleaning $$i/" ; \
	  $(MAKE) clean; \
	  cd .. ; \
	done
