## Process this file with automake to produce Makefile.in

SUBDIRS = src

install-data-local:
	@$(NORMAL_INSTALL)
	if test -d $(srcdir)/pixmaps; then \
	  $(mkinstalldirs) $(DESTDIR)$(pkgdatadir)/pixmaps; \
	  for pixmap in $(srcdir)/pixmaps/*; do \
	    if test -f $$pixmap; then \
	      $(INSTALL_DATA) $$pixmap $(DESTDIR)$(pkgdatadir)/pixmaps; \
	    fi \
	  done \
	fi
	$(INSTALL_DATA) src/gui/siril3.glade $(DESTDIR)$(pkgdatadir)
	$(INSTALL_DATA) src/gui/gtk.css $(DESTDIR)$(pkgdatadir)
	$(mkinstalldirs) $(DESTDIR)$(docdir)
	$(INSTALL_DATA) AUTHORS LICENSE $(DESTDIR)$(docdir)
	$(mkinstalldirs) $(DESTDIR)$(mandir)/man1
	$(INSTALL_DATA) siril.1.gz $(DESTDIR)$(mandir)/man1
