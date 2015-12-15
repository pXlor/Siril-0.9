#ifndef SIRIL_OPENCV_H_
#define SIRIL_OPENCV_H_

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#ifdef HAVE_OPENCV

#ifdef __cplusplus
extern "C" {
#endif

#include "registration/matching/misc.h"

int cvResizeGaussian(fits *, int, int, int);
int cvRotateImage(fits *, double, int, int);
int cvTransformImage(fits *, TRANS, int);
int cvUnsharpFilter(fits*, double, double);
int cvComputeFinestScale(fits *image);
#ifdef __cplusplus
}
#endif

#endif	/* HAVE_OPENCV */

#endif
