#ifndef SIRIL_OPENCV_H_
#define SIRIL_OPENCV_H_

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#ifdef HAVE_OPENCV

#ifdef __cplusplus
extern "C" {
#endif

#include "misc.h"

int resize_gaussian(fits *, int, int, int);
int rotate_image(fits *, double, int, int);
int transforme_image(fits *, TRANS, int);
int unsharp_filter(fits*, double, double);
#ifdef __cplusplus
}
#endif

#endif	/* HAVE_OPENCV */

#endif
