/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2015 team free-astro (see more in AUTHORS file)
 * Reference site is http://free-astro.vinvin.tf/index.php/Siril
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 *
 * Useful links about OpenCV:
 * http://docs.opencv.org/modules/core/doc/intro.html
 * http://docs.opencv.org/modules/imgproc/doc/geometric_transformations.html#resize
 */
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#ifdef HAVE_OPENCV
#include <math.h>
#include <assert.h>
#include "siril.h"
#include "proto.h"
#include "misc.h"
#include "opencv.h"
#include <iostream>
#include <iomanip>
#include "opencv.h"
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace cv;

WORD *fits_to_bgrbgr(fits *image) {
	int ndata = image->rx * image->ry;
	WORD *bgrbgr = new WORD[ndata * 3];
	for (int i = 0, j = 0; i < ndata * 3; i += 3, j++) {
		bgrbgr[i + 0] = image->pdata[BLAYER][j];
		bgrbgr[i + 1] = image->pdata[GLAYER][j];
		bgrbgr[i + 2] = image->pdata[RLAYER][j];
	}
	return bgrbgr;
}

/* resizes image to the sizes toX * toY, and stores it back in image */
int resize_gaussian(fits *image, int toX, int toY, int interpolation) {
	assert(image->data);
	assert(image->rx);

	WORD *bgrbgr = fits_to_bgrbgr(image);

	Mat in(image->ry, image->rx, CV_16UC3, bgrbgr);
	Mat out(toY, toX, CV_16UC3);

	resize(in, out, out.size(), 0, 0, interpolation);

	image->rx = toX;
	image->naxes[0] = toX;
	image->ry = toY;
	image->naxes[1] = toY;
	WORD *newdata = (WORD*) realloc(image->data,
			toX * toY * sizeof(WORD) * image->naxes[2]);
	if (!newdata) {
		free(newdata);
		return 1;
	}
	image->data = newdata;
	Mat channel[3];
	split(out, channel);

	memcpy(image->data, channel[2].data, toX * toY * sizeof(WORD));
	if (image->naxes[2] == 3) {
		memcpy(image->data + toX * toY, channel[1].data,
				toX * toY * sizeof(WORD));
		memcpy(image->data + toX * toY * 2, channel[0].data,
				toX * toY * sizeof(WORD));
	}

	if (image->naxes[2] == 1) {
		image->pdata[0] = image->data;
		image->pdata[1] = image->data;
		image->pdata[2] = image->data;
	} else {
		image->pdata[0] = image->data;
		image->pdata[1] = image->data + (toX * toY);
		image->pdata[2] = image->data + (toX * toY) * 2;
	}
	delete[] bgrbgr;
	in = Mat();
	out = Mat();
	return 0;
}

/* Rotate an image with the angle "angle" */
int rotate_image(fits *image, double angle, int interpolation, int cropped) {
	assert(image->data);
	assert(image->rx);
	assert(image->ry);

	int ndata = image->rx * image->ry;

	WORD *bgrbgr = fits_to_bgrbgr(image);

	Mat in(image->ry, image->rx, CV_16UC3, bgrbgr);
	Mat out(image->ry, image->rx, CV_16UC3);

	if ((angle == 90.0 || angle == 270.0) && interpolation == -1) {	// fast rotation
		transpose(in, out);
		if (angle == 90.0)
			flip(out, out, 0);
		else
			flip(out, out, 1);
	} else {
		Point2f pt(in.cols / 2.0, in.rows / 2.0);// We take the center of the image. Should we pass this in function parameters ?
		Mat r = getRotationMatrix2D(pt, angle, 1.0);
		if (cropped == 1) {
			warpAffine(in, out, r, in.size(), interpolation);
		} else {

			// determine bounding rectangle
			Rect frame = RotatedRect(pt, in.size(), angle).boundingRect();
			// adjust transformation matrix
			r.at<double>(0, 2) += frame.width / 2.0 - pt.x;
			r.at<double>(1, 2) += frame.height / 2.0 - pt.y;

			warpAffine(in, out, r, frame.size(), interpolation);
			ndata = out.cols * out.rows;
			WORD *newdata = (WORD*) realloc(image->data,
					ndata * image->naxes[2] * sizeof(WORD));
			if (!newdata) {
				free(newdata);
				return 1;
			}
			image->data = newdata;
		}
	}
	Mat channel[3];
	split(out, channel);

	memcpy(image->data, channel[2].data, ndata * sizeof(WORD));
	if (image->naxes[2] == 3) {
		memcpy(image->data + ndata, channel[1].data, ndata * sizeof(WORD));
		memcpy(image->data + ndata * 2, channel[0].data, ndata * sizeof(WORD));
	}

	if (image->naxes[2] == 1) {
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data;
		image->pdata[BLAYER] = image->data;
	} else {
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data + ndata;
		image->pdata[BLAYER] = image->data + ndata * 2;
	}
	image->rx = out.cols;
	image->ry = out.rows;
	image->naxes[0] = image->rx;
	image->naxes[1] = image->ry;

	/* free data */
	delete[] bgrbgr;
	in = Mat();
	out = Mat();
	return 0;
}

int transforme_image(fits *image, TRANS trans, int interpolation) {
	assert(image->data);
	assert(image->rx);
	assert(image->ry);

	int ndata = image->rx * image->ry;

	WORD *bgrbgr = fits_to_bgrbgr(image);

	Mat in(image->ry, image->rx, CV_16UC3, bgrbgr);
	Mat out(image->ry, image->rx, CV_16UC3);

	Point2f pt(0, in.rows);
	double angle = atan(trans.c / trans.b) * 180 / M_PI;
	Mat r = getRotationMatrix2D(pt, -angle, 1.0);
	warpAffine(in, out, r, in.size(), interpolation);

	r = (Mat_<double>(2, 3) << 1, 0, trans.a, 0, 1, -trans.d);
	warpAffine(out, out, r, in.size(), interpolation);

	Mat channel[3];
	split(out, channel);

	memcpy(image->data, channel[2].data, ndata * sizeof(WORD));
	if (image->naxes[2] == 3) {
		memcpy(image->data + ndata, channel[1].data, ndata * sizeof(WORD));
		memcpy(image->data + ndata * 2, channel[0].data, ndata * sizeof(WORD));
	}

	if (image->naxes[2] == 1) {
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data;
		image->pdata[BLAYER] = image->data;
	} else {
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data + ndata;
		image->pdata[BLAYER] = image->data + ndata * 2;
	}
	image->rx = out.cols;
	image->ry = out.rows;
	image->naxes[0] = image->rx;
	image->naxes[1] = image->ry;

	/* free data */
	delete[] bgrbgr;
	in = Mat();
	out = Mat();
	r = Mat();
	return 0;
}

int unsharp_filter(fits* image, double sigma, double amount) {
	assert(image->data);
	assert(image->rx);
	int type = CV_16U;
	if (image->naxes[2] != 1)
		type = CV_16UC3;

	Mat in(image->ry, image->rx, type, image->data);
	Mat out, contrast;
	GaussianBlur(in, out, Size(), sigma);
	if (fabs(amount) > 0.0) {
		Mat sharpened = in * (1 + amount) + out * (-amount);
		out = sharpened.clone();
		sharpened.release();
	}

	memcpy(image->data, out.data,
			image->rx * image->ry * sizeof(WORD) * image->naxes[2]);
	if (image->naxes[2] == 1) {
		image->pdata[0] = image->data;
		image->pdata[1] = image->data;
		image->pdata[2] = image->data;
	} else {
		image->pdata[0] = image->data;
		image->pdata[1] = image->data + (image->rx * image->ry);
		image->pdata[2] = image->data + (image->rx * image->ry) * 2;
	}
	in.release();
	return 0;
}

#endif
