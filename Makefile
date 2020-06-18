SCALELIBDIR = scalelib/src
SCALERMDIR  = scale-rm/src
NET2GDIR    = scale-rm/util/netcdf2grads_h
SNODIR      = scale-rm/util/sno
SCALEGMDIR  = scale-gm/src

all: build

build:
	$(MAKE) -C $(SCALELIBDIR) build
	$(MAKE) -C $(SCALELIBDIR) install
	$(MAKE) build-rm build-net2g build-sno build-gm

install: install-rm install-net2g install-sno install-gm

clean:
	$(MAKE) -C $(SCALERMDIR) clean
	$(MAKE) -C $(NET2GDIR) clean
	$(MAKE) -C $(SNODIR) clean
	$(MAKE) -C $(SCALEGMDIR) clean
	$(MAKE) -C $(SCALELIBDIR) clean

distclean:
	$(MAKE) -C $(SCALERMDIR) distclean
	$(MAKE) -C $(NET2GDIR) distclean
	$(MAKE) -C $(SNODIR) distclean
	$(MAKE) -C $(SCALEGMDIR) distclean
	$(MAKE) -C $(SCALELIBDIR) distclean

allclean:
	$(MAKE) -C $(SCALERMDIR) allclean
	$(MAKE) -C $(NET2GDIR) distclean
	$(MAKE) -C $(SNODIR) distclean
	$(MAKE) -C $(SCALEGMDIR) allclean
	$(MAKE) -C $(SCALELIBDIR) allclean

build-rm:
	$(MAKE) -C $(SCALERMDIR) build

build-net2g:
	$(MAKE) -C $(NET2GDIR) build

build-sno:
	$(MAKE) -C $(SNODIR) build

build-gm:
	$(MAKE) -C $(SCALEGMDIR) build

install-rm:
	$(MAKE) -C $(SCALERMDIR) install

install-net2g:
	$(MAKE) -C $(NET2GDIR) install

install-sno:
	$(MAKE) -C $(SNODIR) install

install-gm:
	$(MAKE) -C $(SCALEGMDIR) install
