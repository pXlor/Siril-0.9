## Process this file with automake to produce Makefile.in

bin_PROGRAMS = siril

AM_CPPFLAGS = -DPACKAGE_DATA_DIR='"$(datadir)/${PACKAGE}"' -DPACKAGE_DOC_DIR='"${datarootdir}/doc/${PACKAGE}"' $(GTK_CFLAGS) $(LIBTIFF_CFLAGS) $(LIBPNG_CFLAGS) $(FFTW_CFLAGS) $(CFITSIO_CFLAGS) $(FFMPEG_CFLAGS) $(GSL_CFLAGS) $(FFMS2_CFLAGS) $(LIBCONFIGXX_CFLAGS)

siril_SOURCES = \
	main.c \
	core/siril.c core/siril.h core/command.c core/command.h core/proto.h core/undo.c core/undo.h core/utils.c \
	io/conversion.c io/conversion.h io/ser.c io/ser.h io/films.c io/films.h \
	io/image_formats_libraries.c io/image_formats_internal.c io/image_format_fits.c \
	io/sequence.c io/seqfile.c io/single_image.c io/single_image.h \
	gui/callbacks.c gui/callbacks.h gui/histogram.c gui/histogram.h \
	gui/PSF_list.h gui/PSF_list.c gui/sequence_list.c gui/registration_preview.c\
	registration/registration.c registration/registration.h \
	registration/matching/match.c registration/matching/atpmatch.c registration/matching/misc.c \
	stacking/stacking.c stacking/stacking.h \
	algos/gradient.c algos/gradient.h algos/quality.c algos/quality.h \
	algos/fft.c algos/fft.h \
	algos/colors.c algos/colors.h \
	algos/demosaicing.c algos/demosaicing.h \
	algos/pave.c algos/transform.c algos/io_wave.c algos/reconstr.c \
	algos/Def_Math.h algos/Def_Wavelet.h algos/Def_Mem.h \
	algos/PSF.c algos/PSF.h algos/star_finder.c algos/star_finder.h \
	compositing/compositing.c compositing/compositing.h compositing/filters.c compositing/filters.h
	
siril_LDADD = $(GTK_LIBS) $(LIBRAW_LIBS) $(LIBTIFF_LIBS) $(LIBPNG_LIBS) $(FFTW_LIBS) $(CFITSIO_LIBS) $(GSL_LIBS) $(FFMS2_LIBS) $(LIBCONFIGXX_LIBS)

if HAVE_OPENCV
siril_SOURCES += opencv/opencv.cpp opencv/opencv.h opencv/ecc/ecc.cpp opencv/ecc/ecc.h
siril_LDADD += $(OPENCV_LIBS)
siril_LINK = $(CXXLINK)
else
siril_LINK = $(LINK)
endif
