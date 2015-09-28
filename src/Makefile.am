## Process this file with automake to produce Makefile.in

bin_PROGRAMS = siril

AM_CPPFLAGS = -DPACKAGE_DATA_DIR='"$(datadir)/${PACKAGE}"' -DPACKAGE_DOC_DIR='"${datarootdir}/doc/${PACKAGE}"' $(GTK_CFLAGS) $(LIBTIFF_CFLAGS) $(LIBPNG_CFLAGS) $(FFTW_CFLAGS) $(CFITSIO_CFLAGS) $(FFMPEG_CFLAGS) $(GSL_CFLAGS) $(FFMS2_CFLAGS) $(LIBCONFIGXX_CFLAGS)

siril_SOURCES = \
	siril.c main.c command.c command.h \
	utils.c siril.h proto.h \
	callbacks.c callbacks.h \
	registration.c registration.h registration_preview.c \
	quality.c quality.h \
	stacking.c stacking.h \
	fft.c fft.h\
	colors.c colors.h \
	demosaicing.c demosaicing.h \
	sequence.c seqfile.c sequence_list.c \
	conversion.c \
	single_image.c single_image.h \
	histogram.c histogram.h \
	ser.c ser.h films.c films.h \
	image_formats_libraries.c image_formats_internal.c image_format_fits.c \
	pave.c transform.c io_wave.c reconstr.c Def_Math.h Def_Wavelet.h Def_Mem.h \
	PSF.c PSF.h PSF_list.c PSF_list.h star_finder.c star_finder.h gradient.c gradient.h \
	compositing.c compositing.h filters.c filters.h \
	match.c atpmatch.c misc.c match.h atpmatch.h misc.h \
	undo.c undo.h

siril_LDADD = $(GTK_LIBS) $(LIBRAW_LIBS) $(LIBTIFF_LIBS) $(LIBPNG_LIBS) $(FFTW_LIBS) $(CFITSIO_LIBS) $(GSL_LIBS) $(FFMS2_LIBS) $(LIBCONFIGXX_LIBS)

if HAVE_OPENCV
siril_SOURCES += opencv.cpp opencv.h 
siril_LDADD += $(OPENCV_LIBS)
siril_LINK = $(CXXLINK)
else
siril_LINK = $(LINK)
endif
