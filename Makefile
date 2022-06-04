SCALELIBDIR = scalelib/src
SCALERMDIR  = scale-rm/src
SNODIR      = scale-rm/util/sno
SCALEGMDIR  = scale-gm/src

all: build

build:
	$(MAKE) -C $(SCALELIBDIR) build
	$(MAKE) -C $(SCALELIBDIR) install
	$(MAKE) build-rm build-sno build-gm

install: install-rm install-sno install-gm

clean:
	$(MAKE) -C $(SCALERMDIR) clean
	$(MAKE) -C $(SNODIR) clean
	$(MAKE) -C $(SCALEGMDIR) clean
	$(MAKE) -C $(SCALELIBDIR) clean

distclean:
	$(MAKE) -C $(SCALERMDIR) distclean
	$(MAKE) -C $(SNODIR) distclean
	$(MAKE) -C $(SCALEGMDIR) distclean
	$(MAKE) -C $(SCALELIBDIR) distclean

allclean:
	$(MAKE) -C $(SCALERMDIR) allclean
	$(MAKE) -C $(SNODIR) distclean
	$(MAKE) -C $(SCALEGMDIR) allclean
	$(MAKE) -C $(SCALELIBDIR) allclean

build-rm:
	$(MAKE) -C $(SCALERMDIR) build

build-sno:
	$(MAKE) -C $(SNODIR) build

build-gm:
	$(MAKE) -C $(SCALEGMDIR) build

install-rm:
	$(MAKE) -C $(SCALERMDIR) install

install-sno:
	$(MAKE) -C $(SNODIR) install

install-gm:
	$(MAKE) -C $(SCALEGMDIR) install
