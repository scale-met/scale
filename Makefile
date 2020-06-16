INSTALL=install

all: build-rm build-net2g

install: install-rm install-net2g

clean:
	cd scale-rm/src && $(MAKE) clean allclean
	cd scale-rm/util/netcdf2grads_h && $(MAKE) clean
	cd scale-gm/src && $(MAKE) clean allclean

build-rm:
	cd scale-rm/src && $(MAKE)

build-net2g:
	cd scale-rm/util/netcdf2grads_h && $(MAKE)

build-gm:
	cd scale-gm/src && $(MAKE)

bin/scale-rm : build-rm
bin/net2g : build-net2g
bin/scale-gm : build-gm

install-rm: bin/scale-rm
	$(INSTALL) -d $(PREFIX)/bin
	$(INSTALL) bin/scale-rm $(PREFIX)/bin
	$(INSTALL) bin/scale-rm_init $(PREFIX)/bin
	$(INSTALL) bin/scale-rm_pp $(PREFIX)/bin

install-net2g: bin/net2g
	$(INSTALL) bin/net2g $(PREFIX)/bin

install-gm: bin/scale-gm
	$(INSTALL) -d $(PREFIX)/bin
	$(INSTALL) bin/scale-gm $(PREFIX)/bin
	$(INSTALL) bin/gm_fio_cat $(PREFIX)/bin
	$(INSTALL) bin/gm_fio_dump $(PREFIX)/bin
	$(INSTALL) bin/gm_fio_ico2ll $(PREFIX)/bin
	$(INSTALL) bin/gm_fio_sel $(PREFIX)/bin
	$(INSTALL) bin/gm_mkhgrid $(PREFIX)/bin
	$(INSTALL) bin/gm_mkllmap $(PREFIX)/bin
	$(INSTALL) bin/gm_mkmnginfo $(PREFIX)/bin
	$(INSTALL) bin/gm_mkrawgrid $(PREFIX)/bin
	$(INSTALL) bin/gm_mkvlayer $(PREFIX)/bin
