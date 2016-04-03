dnl Process this file with autoconf to produce a configure script.

AC_INIT(siril, 0.9.2)
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_SRCDIR([src/main.c])

remember_set_CFLAGS="$CFLAGS"

AC_PROG_CC
AC_HEADER_STDC
AC_OPENMP

dnl remove -g flag from default CFLAGS
if test "x$remember_set_CFLAGS" = "x"; then
	if test "$CFLAGS" = "-g -O2"; then
		CFLAGS="-O2"
	elif test "$CFLAGS" = "-g"; then
		CFLAGS=""
	fi
fi

CFLAGS="$CFLAGS $OPENMP_CFLAGS"
CXXFLAGS="$CFLAGS"

if test `uname` != "Darwin"; then
	LDFLAGS="$LDFLAGS -rdynamic"
fi

dnl Set PACKAGE_DATA_DIR in config.h.
dnl if test "x${datadir}" = 'x${prefix}/share'; then
dnl   if test "x${prefix}" = "xNONE"; then
dnl     AC_DEFINE_UNQUOTED(PACKAGE_DATA_DIR, "${ac_default_prefix}/share/${PACKAGE}", [Package data dir1])
dnl   else
dnl     AC_DEFINE_UNQUOTED(PACKAGE_DATA_DIR, "${prefix}/share/${PACKAGE}", [Package data dir2])
dnl   fi
dnl else
dnl   AC_DEFINE_UNQUOTED(PACKAGE_DATA_DIR, "${datadir}/${PACKAGE}", [Package data dir3])
dnl fi

dnl Set PACKAGE_SOURCE_DIR in config.h.
packagesrcdir=`cd $srcdir && pwd`
AC_DEFINE_UNQUOTED(PACKAGE_SOURCE_DIR, "${packagesrcdir}", [Package source dir])

dnl Set PACKAGE_DOC_DIR in config.h.
dnl AC_DEFINE_UNQUOTED(PACKAGE_DOC_DIR, "${ac_default_prefix}/share/doc/${PACKAGE}", "package doc dir")

dnl ================ MANDATORY DEPENDECIES ================ 
dnl check math lib
AC_CHECK_LIB(m, log, LIBS="$LIBS -lm")

dnl check GTK
PKG_CHECK_MODULES(GTK, gtk+-3.0 >= 3.6.0)

dnl check fftw3
PKG_CHECK_MODULES(FFTW, [fftw3])

dnl check GNU Scientific Library
PKG_CHECK_MODULES(GSL, [gsl < 2],
		  [AC_DEFINE([HAVE_GSL_1], [1], [Use GSL 1])],
		  [PKG_CHECK_MODULES(GSL, [gsl >= 2],
				     [AC_DEFINE([HAVE_GSL_2], [1], [Use GSL 2])
				     ])
		  ])

dnl check cfitsio
PKG_CHECK_MODULES(CFITSIO, [cfitsio])

dnl check libconfig
PKG_CHECK_MODULES([LIBCONFIGXX], [libconfig++ >= 1.4],,
	AC_MSG_ERROR([libconfig++ 1.4 or newer not found.]))

dnl ================ OPTIONAL DEPENDECIES ================ 
dnl check raw lib
PKG_CHECK_MODULES(LIBRAW, [libraw],
 		  [AC_DEFINE([HAVE_LIBRAW], [1], [Using RAW images])],
 		  AC_MSG_WARN([libraw not found. Not using RAW importer.]))

dnl check tiff lib
PKG_CHECK_MODULES(LIBTIFF, [libtiff-4],
		  [AC_DEFINE([HAVE_LIBTIFF], [1], [Using TIFF images])],
		  AC_MSG_WARN([libtiff not found. Not using TIFF importer and exporter.]))

dnl check jpeg lib
dnl PKG_CHECK_MODULES(JPEG, [libjpeg])
AC_CHECK_LIB(jpeg, jpeg_mem_src, [],
	     AC_MSG_WARN([libjpeg not found. Not using JPEG importer and exporter.]))

dnl check png lib
PKG_CHECK_MODULES(LIBPNG, [libpng],
		  [AC_DEFINE([HAVE_LIBPNG], [1], [Using PNG images])],
		  AC_MSG_WARN([libpng not found. Not using PNG importer and exporter.]))
		  
dnl check gif lib
AC_CHECK_LIB(gif, DGifOpenFileHandle, [])
AS_IF([test "$ac_cv_lib_gif_DGifOpenFileHandle" != yes],
      AC_MSG_WARN([giflib not found. Not using GIF exporter.]))

dnl check ffms2		  
PKG_CHECK_MODULES(FFMS2, [ffms2 < 2.21],
		  [AC_DEFINE([HAVE_FFMS2_1], [1], [Use FFMS2 < 2.21])],
		  [PKG_CHECK_MODULES(FFMS2, [ffms2 >= 2.21],
				     [AC_DEFINE([HAVE_FFMS2_2], [1], [Use FFMS2 > 2.21])
				     ])
		  ])

dnl check opencv
AC_PROG_CXX
if test -n "$CXX" && ( test "X$CXX" != "Xno" &&
		      ( (test "X$CXX" = "Xg++" && `g++ -v >/dev/null 2>&1` ) ||
		      (test "X$CXX" != "Xg++"))) ; then
	AC_PROG_CXXCPP
	dnl AC_MSG_WARN(CXX is "$CXX")
	dnl http://nerdland.net/2009/07/detecting-c-libraries-with-autotools/ solution 2
	PKG_CHECK_MODULES([OPENCV], [opencv >= 2.0],
			  [AC_DEFINE([HAVE_OPENCV], [1], [Using OpenCV])]
			  [OPENCV_LIBS="-lopencv_core -lopencv_imgproc"],
			  AC_MSG_WARN([opencv not found. Not using some image processing.]))
else
	AC_MSG_WARN([no C++ compiler found, opencv cannot be used])
fi
AM_CONDITIONAL(HAVE_OPENCV, [test "x$OPENCV_LIBS" != "x"])

dnl if test "x$CXX" = "xno"; then
dnl 	dnl TODO: is it possible to require CXX only when there is opencv?
dnl 	AC_MSG_WARN([no C++ compiler found, opencv cannot be used])
dnl else
dnl 	AC_MSG_WARN(CXX is "$CXX")
dnl 	dnl http://nerdland.net/2009/07/detecting-c-libraries-with-autotools/ solution 2
dnl PKG_CHECK_MODULES([OPENCV], [opencv >= 2.0],
dnl  		  [AC_DEFINE([HAVE_OPENCV], [1], [Using OpenCV])]
dnl 		  [OPENCV_LIBS="-lopencv_core -lopencv_imgproc"],
dnl 		  AC_MSG_WARN([opencv not found. Not using image resizing.]))
dnl fi

dnl checking several ffmpeg libraries - replaced by ffms2
dnl can replace ffms2, but cannot have both enabled at the same time
dnl PKG_CHECK_MODULES(FFMPEG, [libavformat libavutil libavcodec libswscale],
dnl 		  [AC_DEFINE([HAVE_FFMPEG], [1], [Using AVI files natively])],
dnl 		  AC_MSG_WARN([ffmpeg libav libraries not found. Not using AVI format.]))

dnl netpbm is not used anymore, it is replaced by an internal PPM/PNM loader
dnl AC_CHECK_LIB(netpbm, pgm_init, netpbm="true",
dnl   AC_MSG_WARN([libnetpbm not found. jpeg and png conversions will not be available.]))
dnl AM_CONDITIONAL(HAVE_NETPBM, test x$netpbm = xtrue)
dnl ====================================================== 

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

