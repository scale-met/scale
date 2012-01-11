SUBDIRS=boundary grid initial refstate

all:
	@for i in $(SUBDIRS) ; do \
	  if [ -d $$i ]; then \
	      (cd $$i ; \
	      echo "Making $$i/" ; \
	      $(MAKE)) ; \
	  fi; \
	done

clean:
	@for i in $(SUBDIRS) ; do \
	  if [ -d $$i ]; then \
	      (cd $$i ; \
	      echo "Cleaning $$i/" ; \
	      $(MAKE) clean); \
	  fi; \
	done
