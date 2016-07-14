/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2016 team free-astro (see more in AUTHORS file)
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
 *
 * WARNING: the code in this file and its header will not work properly
 * on big endian systems.
 */

#define __USE_XOPEN
#define _GNU_SOURCE

#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "core/siril.h"
#include "core/proto.h"
#include "gui/callbacks.h"
#include "algos/demosaicing.h"
#include "io/ser.h"

/* 62135596800 sec from year 0001 to 01 janv. 1970 00:00:00 GMT */
const uint64_t epochTicks = 621355968000000000UL;
const uint64_t ticksPerSecond = 10000000;

int display_date(uint64_t date, char *txt) {
	struct tm *timeinfo;
	char str[256];
	time_t t_sec;

	if (date <= 0)
		return -1;
	t_sec = (time_t) (date - epochTicks) / ticksPerSecond;
	if (t_sec < 0)
		return -1;
	timeinfo = gmtime(&t_sec);
	strcpy(str, txt);
	strftime(str + strlen(txt), 255, "%F %r", timeinfo);
	puts(str);
	return 0;
}

/* Comes from http://linux.die.net/man/3/timegm
 * This is not thread-safe
 */
static time_t __timegm(struct tm *tm) {
	time_t ret;
	char *tz;

	tz = getenv("TZ");
	setenv("TZ", "", 1);
	tzset();
	ret = mktime(tm);
	if (tz)
		setenv("TZ", tz, 1);
	else
		unsetenv("TZ");
	tzset();
	return ret;
}

int convert_char_to_time(char *date, uint64_t *utc, uint64_t *local) {
	struct tm timeinfo = { };
	time_t ut, t;

	if (date[0] == '\0')
		return -1;
	strptime(date, "%FT%T", &timeinfo); /* YYYY-MM-DDThh:mm:ss */
	timeinfo.tm_isdst = -1;

	/* get UTC time from timeinfo* */
	ut = __timegm(&timeinfo);
	ut *= ticksPerSecond;
	ut += epochTicks;
	*utc = (uint64_t) ut;

	/* get local time from timeinfo* */
	t = mktime(&timeinfo);
	t *= ticksPerSecond;
	t += epochTicks;
	*local = (uint64_t) t;

	return 0;
}

char *convert_color_id_to_char(ser_color color_id) {
	switch (color_id) {
	case SER_MONO:
		return "MONO";
	case SER_BAYER_RGGB:
		return "RGGB";
	case SER_BAYER_BGGR:
		return "BGGR";
	case SER_BAYER_GBRG:
		return "GBRG";
	case SER_BAYER_GRBG:
		return "GRBG";
	case SER_BAYER_CYYM:
		return "CYYM";
	case SER_BAYER_YCMY:
		return "YCMY";
	case SER_BAYER_YMCY:
		return "YMCY";
	case SER_BAYER_MYYC:
		return "MYYC";
	case SER_RGB:
		return "RGB";
	case SER_BGR:
		return "BGR";
	default:
		return "";
	}
}

void ser_display_info(struct ser_struct *ser_file) {
	char *color = convert_color_id_to_char(ser_file->color_id);

	fprintf(stdout, "=========== SER file info ==============\n");
	fprintf(stdout, "file id: %s\n", ser_file->file_id);
	fprintf(stdout, "lu id: %d\n", ser_file->lu_id);
	fprintf(stdout, "sensor type: %s\n", color);
	fprintf(stdout, "image size: %d x %d (%d bits)\n", ser_file->image_width,
			ser_file->image_height, ser_file->bit_pixel_depth);
	fprintf(stdout, "frame count: %u\n", ser_file->frame_count);
	fprintf(stdout, "observer: %s\n", ser_file->observer);
	fprintf(stdout, "instrument: %s\n", ser_file->instrument);
	fprintf(stdout, "telescope: %s\n", ser_file->telescope);
	display_date(ser_file->date, "local time: ");
	display_date(ser_file->date_utc, "UTC time: ");
	fprintf(stdout, "========================================\n");
}

int _ser_read_header(struct ser_struct *ser_file) {
	char header[SER_HEADER_LEN];
	if (!ser_file || ser_file->fd <= 0)
		return -1;
	if ((off_t) -1 == lseek(ser_file->fd, 0, SEEK_SET)) {
		perror("seek");
		return -1;
	}
	if (sizeof(header) != read(ser_file->fd, header, sizeof(header))) {
		perror("read");
		return -1;
	}
	// modify this to support big endian
	memcpy(&ser_file->lu_id, header + 14, 28);	// read all integers
	memcpy(&ser_file->date, header + 162, 16);

	// strings
	ser_file->file_id = strndup(header, 14);
	ser_file->observer = strndup(header + 42, 40);
	ser_file->instrument = strndup(header + 82, 40);
	ser_file->telescope = strndup(header + 122, 40);

	/* internal representations of header data */
	if (ser_file->bit_pixel_depth <= 8)
		ser_file->byte_pixel_depth = SER_PIXEL_DEPTH_8;
	else ser_file->byte_pixel_depth = SER_PIXEL_DEPTH_16;

	if (ser_file->color_id == SER_RGB || ser_file->color_id == SER_BGR)
		ser_file->number_of_planes = 3;
	else ser_file->number_of_planes = 1;

	return 0;
}

int ser_write_header(struct ser_struct *ser_file) {
	char header[SER_HEADER_LEN];
	if (!ser_file || ser_file->fd <= 0)
		return -1;
	if ((off_t) -1 == lseek(ser_file->fd, 0, SEEK_SET)) {
		perror("seek");
		return -1;
	}

	// modify this to support big endian
	memset(header, 0, sizeof(header));
	memcpy(header, ser_file->file_id, 14);
	memcpy(header + 14, &ser_file->lu_id, 28);
	memcpy(header + 42, ser_file->observer, 40);
	memcpy(header + 82, ser_file->instrument, 40);
	memcpy(header + 122, ser_file->telescope, 40);
	memcpy(header + 162, &ser_file->date, 16);

	if (sizeof(header) != write(ser_file->fd, header, sizeof(header))) {
		perror("write");
		return 1;
	}
	return 0;
}

int ser_write_and_close(struct ser_struct *ser_file) {
	ser_write_header(ser_file);	// writes the header
	return ser_close_file(ser_file);// closes, frees and zeroes
}

/* ser_file must be allocated */
int ser_create_file(const char *filename, struct ser_struct *ser_file, gboolean overwrite, struct ser_struct *copy_from) {
	if (overwrite)
		unlink(filename);
	if ((ser_file->fd = open(filename, O_CREAT | O_RDWR,
			S_IWRITE | S_IREAD)) == -1) {
		perror("open SER file for creation");
		return 1;
	}

	if (copy_from) {
		memcpy(&ser_file->lu_id, &copy_from->lu_id, 28);
		memcpy(&ser_file->date, &copy_from->date, 16);
		memcpy(&ser_file->date_utc, &copy_from->date_utc, 16);
		ser_file->file_id = strdup(copy_from->file_id);
		ser_file->observer = strdup(copy_from->observer);
		ser_file->instrument = strdup(copy_from->instrument);
		ser_file->telescope = strdup(copy_from->telescope);
		ser_file->byte_pixel_depth = copy_from->byte_pixel_depth;
		ser_file->number_of_planes = copy_from->number_of_planes;
		/* we write the header now, but it should be written again
		 * before closing in case the number of the image in the new
		 * SER changes from the copied SER */
		ser_write_header(ser_file);
	} else {	// new SER
		ser_file->file_id = strdup("Made by Siril");
		ser_file->lu_id = 0;
		ser_file->little_endian = SER_LITTLE_ENDIAN; // what will it do on big endian machine?
		ser_file->observer = strdup("");
		ser_file->instrument = strdup("");
		ser_file->telescope = strdup("");
		memset(&ser_file->date, 0, 16);
		ser_file->number_of_planes = 0;	// used as an indicator of new SER

		/* next operation should be ser_write_frame_from_fit, which writes with no
		 * seek and expects to be after the header */
		if ((off_t) -1 == lseek(ser_file->fd, SER_HEADER_LEN, SEEK_SET)) {
			perror("seek");
			return -1;
		}
	}
	ser_file->filename = strdup(filename);
	ser_file->frame_count = 0;	// incremented on image add
#ifdef _OPENMP
	omp_init_lock(&ser_file->fd_lock);
#endif
	siril_log_message("Created SER file %s\n", filename);
	return 0;
}

/* populate fields that are not already set in ser_create_file */
void ser_header_from_fit(struct ser_struct *ser_file, fits *fit) {
	ser_file->image_width = fit->rx;
	ser_file->image_height = fit->ry;
	fprintf(stdout, "setting SER image size as %dx%d\n", fit->rx, fit->ry);
	if (fit->naxes[2] == 1) {
		// maybe it's read as CFA...
		ser_file->color_id = SER_MONO;
	} else if (fit->naxes[2] == 3) {
		ser_file->color_id = SER_RGB;
	}
	if (ser_file->color_id == SER_RGB || ser_file->color_id == SER_BGR)
		ser_file->number_of_planes = 3;
	else ser_file->number_of_planes = 1;

	if (fit->bitpix == BYTE_IMG) {
		ser_file->byte_pixel_depth = SER_PIXEL_DEPTH_8;
		ser_file->bit_pixel_depth = 8;
	} else if (fit->bitpix == USHORT_IMG || fit->bitpix == SHORT_IMG) {
		ser_file->byte_pixel_depth = SER_PIXEL_DEPTH_16;
		ser_file->bit_pixel_depth = 16;
	} else {
		siril_log_message("Writing to SER files from larger than 16-bit FITS images is not yet implemented\n");
	}
	if (fit->instrume[0] != 0) {
		free(ser_file->instrument);
		ser_file->instrument = strdup(fit->instrume);
	}
	int ret = convert_char_to_time(fit->date_obs, &ser_file->date_utc, &ser_file->date);
	if (ret == -1)
		convert_char_to_time(fit->date, &ser_file->date_utc, &ser_file->date);
	// TODO: copy data from the fit header: observer, telescope
}

int ser_open_file(char *filename, struct ser_struct *ser_file) {
	if (ser_file->fd > 0) {
		fprintf(stderr, "SER: file already opened, or badly closed\n");
		return -1;
	}
	ser_file->fd = open(filename, O_RDONLY);
	if (ser_file->fd == -1) {
		perror("SER file open");
		return -1;
	}
	if (_ser_read_header(ser_file)) {
		fprintf(stderr, "SER: reading header failed, closing file %s\n",
				filename);
		ser_close_file(ser_file);
		return -1;
	}
	ser_file->filename = strdup(filename);

#ifdef _OPENMP
	omp_init_lock(&ser_file->fd_lock);
#endif
	return 0;
}

int ser_close_file(struct ser_struct *ser_file) {
	int retval = 0;
	if (!ser_file)
		return retval;
	if (ser_file->fd > 0) {
		retval = close(ser_file->fd);
		ser_file->fd = -1;
		// removing the created file if frames have not been added
		if (ser_file->frame_count < 2 && ser_file->filename) {
			unlink(ser_file->filename);
			retval = 1;
		}
	}
	if (ser_file->file_id)
		free(ser_file->file_id);
	if (ser_file->observer)
		free(ser_file->observer);
	if (ser_file->instrument)
		free(ser_file->instrument);
	if (ser_file->telescope)
		free(ser_file->telescope);
	if (ser_file->filename)
		free(ser_file->filename);
#ifdef _OPENMP
	omp_destroy_lock(&ser_file->fd_lock);
#endif
	ser_init_struct(ser_file);
	return retval;
}

void ser_init_struct(struct ser_struct *ser_file) {
	memset(ser_file, 0, sizeof(struct ser_struct));
}

/* frame number starts at 0 */
int ser_read_frame(struct ser_struct *ser_file, int frame_no, fits *fit) {
	int retval, frame_size, i, j, swap = 0;
	off_t offset;
	WORD *olddata, *tmp;
	if (!ser_file || ser_file->fd <= 0 || !ser_file->number_of_planes ||
			!fit || frame_no < 0 || frame_no >= ser_file->frame_count)
		return -1;
	frame_size = ser_file->image_width * ser_file->image_height
			* ser_file->number_of_planes;
	olddata = fit->data;
	if ((fit->data = realloc(fit->data, frame_size * sizeof(WORD))) == NULL) {
		fprintf(stderr, "ser_read: error realloc %s %d\n", ser_file->filename,
				frame_size);
		if (olddata)
			free(olddata);
		return -1;
	}

	offset = SER_HEADER_LEN + (off_t)frame_size *
		(off_t)ser_file->byte_pixel_depth * (off_t)frame_no;
	/*fprintf(stdout, "offset is %lu (frame %d, %d pixels, %d-byte)\n", offset,
	 frame_no, frame_size, ser_file->pixel_bytedepth);*/
#ifdef _OPENMP
	omp_set_lock(&ser_file->fd_lock);
#endif
	if ((off_t) -1 == lseek(ser_file->fd, offset, SEEK_SET)) {
#ifdef _OPENMP
		omp_unset_lock(&ser_file->fd_lock);
#endif
		return -1;
	}
	retval = read(ser_file->fd, fit->data, frame_size * ser_file->byte_pixel_depth);
#ifdef _OPENMP
	omp_unset_lock(&ser_file->fd_lock);
#endif
	if (retval != frame_size * ser_file->byte_pixel_depth)
		return -1;

	ser_manage_endianess_and_depth(ser_file, fit->data, frame_size);

	fit->bitpix = (ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_8) ? BYTE_IMG : USHORT_IMG;

	/* If the user checks the SER CFA box, the video is opened in B&W
	 * RGB and BGR are not coming from raw data. In consequence CFA does
	 * not exist for these kind of cam */
	ser_color type_ser = ser_file->color_id;
	if (!com.debayer.open_debayer && type_ser != SER_RGB && type_ser != SER_BGR)
		type_ser = SER_MONO;

	switch (type_ser) {
	case SER_MONO:
		fit->naxis = 2;
		fit->naxes[0] = fit->rx = ser_file->image_width;
		fit->naxes[1] = fit->ry = ser_file->image_height;
		fit->naxes[2] = 1;
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->data;
		fit->pdata[BLAYER] = fit->data;
		break;
	case SER_BAYER_RGGB:
	case SER_BAYER_BGGR:
	case SER_BAYER_GBRG:
	case SER_BAYER_GRBG:
		fit->naxes[0] = fit->rx = ser_file->image_width;
		fit->naxes[1] = fit->ry = ser_file->image_height;
		if (com.debayer.ser_use_bayer_header)
			set_combo_box_bayer_pattern(type_ser);
		debayer(fit, com.debayer.bayer_inter);
		break;
	case SER_BGR:
		swap = 2;
		/* no break */
	case SER_RGB:
		tmp = malloc(frame_size * sizeof(WORD));
		memcpy(tmp, fit->data, sizeof(WORD) * frame_size);
		fit->naxes[0] = fit->rx = ser_file->image_width;
		fit->naxes[1] = fit->ry = ser_file->image_height;
		fit->naxes[2] = 3;
		fit->naxis = 3;
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->data + fit->rx * fit->ry;
		fit->pdata[BLAYER] = fit->data + fit->rx * fit->ry * 2;
		for (i = 0, j = 0; j < fit->rx * fit->ry; i += 3, j++) {
			fit->pdata[0 + swap][j] = tmp[i + RLAYER];
			fit->pdata[1       ][j] = tmp[i + GLAYER];
			fit->pdata[2 - swap][j] = tmp[i + BLAYER];
		}
		free(tmp);
		break;
	case SER_BAYER_CYYM:
	case SER_BAYER_YCMY:
	case SER_BAYER_YMCY:
	case SER_BAYER_MYYC:
	default:
		siril_log_message("This type of Bayer pattern is not handled yet.\n");
		return -1;
	}
	fits_flip_top_to_bottom(fit);
	return 0;
}

/* read an area of an image in an opened SER sequence */
int ser_read_opened_partial(struct ser_struct *ser_file, int layer,
		int frame_no, WORD *buffer, const rectangle *area) {
	off_t offset;
	int frame_size, read_size, retval, xoffset, yoffset, x, y, color_offset;
	ser_color type_ser;
	WORD *rawbuf, *demosaiced_buf, *rgb_buf;
	rectangle debayer_area, image_area;

	if (!ser_file || ser_file->fd <= 0 || frame_no < 0
			|| frame_no >= ser_file->frame_count)
		return -1;
	frame_size = ser_file->image_width * ser_file->image_height *
		ser_file->number_of_planes * ser_file->byte_pixel_depth;

	/* If the user checks the SER CFA box, the video is opened in B&W
	 * RGB and BGR are not coming from raw data. In consequence CFA does
	 * not exist for these kind of cam */
	type_ser = ser_file->color_id;
	if (!com.debayer.open_debayer && type_ser != SER_RGB && type_ser !=
			SER_BGR)
		type_ser = SER_MONO;

	switch (type_ser) {
	case SER_MONO:
		offset = SER_HEADER_LEN + (off_t)frame_size * (off_t)frame_no +	// requested frame
			(off_t)(area->y * ser_file->image_width + area->x)
						* ser_file->byte_pixel_depth;	// requested area
#ifdef _OPENMP
		omp_set_lock(&ser_file->fd_lock);
#endif
		if ((off_t) -1 == lseek(ser_file->fd, offset, SEEK_SET)) {
#ifdef _OPENMP
			omp_unset_lock(&ser_file->fd_lock);
#endif
			return -1;
		}
		read_size = area->w * area->h * ser_file->byte_pixel_depth;
		retval = read(ser_file->fd, buffer, read_size);
#ifdef _OPENMP
		omp_unset_lock(&ser_file->fd_lock);
#endif
		if (retval != read_size)
			return -1;

		ser_manage_endianess_and_depth(ser_file, buffer, area->w * area->h);
		break;

	case SER_BAYER_RGGB:
	case SER_BAYER_BGGR:
	case SER_BAYER_GBRG:
	case SER_BAYER_GRBG:
		/* SER v2: RGB images obtained from demosaicing.
		 * Original is monochrome, we demosaic it in an area slightly larger than the
		 * requested area, giving 3 channels in form of RGBRGBRGB buffers, and finally
		 * we extract one of the three channels and crop it to the requested area. */
		if (com.debayer.ser_use_bayer_header)
			set_combo_box_bayer_pattern(type_ser);
		if (layer < 0 || layer >= 3) {
			siril_log_message("For a demosaiced image, layer has to be R, G or B (0 to 2).\n");
			return -1;
		}

		image_area = (rectangle) { .x = 0, .y = 0,
			.w = ser_file->image_width, .h = ser_file->image_height };
		get_debayer_area(area, &debayer_area, &image_area, &xoffset, &yoffset);

		offset = SER_HEADER_LEN + frame_size * frame_no +	// requested frame
				(debayer_area.y * ser_file->image_width + debayer_area.x)
						* ser_file->byte_pixel_depth;	// requested area
#ifdef _OPENMP
		omp_set_lock(&ser_file->fd_lock);
#endif
		if ((off_t) -1 == lseek(ser_file->fd, offset, SEEK_SET)) {
#ifdef _OPENMP
			omp_unset_lock(&ser_file->fd_lock);
#endif
			return -1;
		}

		// allocating a buffer for WORD because it's going to be converted in-place
		rawbuf = malloc(debayer_area.w * debayer_area.h * sizeof(WORD));
		if (!rawbuf) {
#ifdef _OPENMP
			omp_unset_lock(&ser_file->fd_lock);
#endif
			siril_log_message("Out of memory - aborting\n");
			return -1;
		}
		read_size = debayer_area.w * debayer_area.h * ser_file->byte_pixel_depth;
		retval = read(ser_file->fd, rawbuf, read_size);
#ifdef _OPENMP
		omp_unset_lock(&ser_file->fd_lock);
#endif
		if (retval != read_size) {
			free(rawbuf);
			perror("read");
			return -1;
		}

		ser_manage_endianess_and_depth(ser_file, rawbuf, debayer_area.w * debayer_area.h);

		demosaiced_buf = debayer_buffer(rawbuf, &debayer_area.w,
				&debayer_area.h, com.debayer.bayer_inter,
				com.debayer.bayer_pattern);
		free(rawbuf);
		if (demosaiced_buf == NULL) {
			return -1;
		}

		/* area is the destination area.
		 * debayer_area is the demosaiced buf area.
		 * xoffset and yoffset are the x,y offsets of area in the debayer area.
		 */
		for (y = 0; y < area->h; y++) {
			for (x = 0; x < area->w; x++) {
				buffer[y*area->w + x] = demosaiced_buf[(yoffset+y)*debayer_area.w*3 + xoffset+x*3 + layer]; 
			}
		}

		free(demosaiced_buf);
		break;
	case SER_BGR:
	case SER_RGB:
		//siril_log_message("Work in progress... Available soon !!\n");
		//return -1;
		assert(ser_file->number_of_planes == 3);

		offset = SER_HEADER_LEN + frame_size * frame_no +	// requested frame
			(area->y * ser_file->image_width + area->x) *
			ser_file->byte_pixel_depth * 3;	// requested area
#ifdef _OPENMP
		omp_set_lock(&ser_file->fd_lock);
#endif
		if ((off_t) -1 == lseek(ser_file->fd, offset, SEEK_SET)) {
#ifdef _OPENMP
			omp_unset_lock(&ser_file->fd_lock);
#endif
			return -1;
		}

		read_size = area->w * area->h * ser_file->byte_pixel_depth * 3;
		// allocating a buffer for WORD because it's going to be converted in-place
		rgb_buf = malloc(area->w * area->h * 3 * sizeof(WORD));
		if (!rgb_buf) {
#ifdef _OPENMP
			omp_unset_lock(&ser_file->fd_lock);
#endif
			siril_log_message("Out of memory - aborting\n");
			return -1;
		}
		retval = read(ser_file->fd, rgb_buf, read_size);
#ifdef _OPENMP
		omp_unset_lock(&ser_file->fd_lock);
#endif
		if (retval != read_size) {
			free(rgb_buf);
			perror("read");
			return -1;
		}

		ser_manage_endianess_and_depth(ser_file, rgb_buf, area->w * area->h * 3);

		color_offset = layer;
		if (type_ser == SER_BGR) {
			color_offset = 2 - layer;
		}

		for (y = 0; y < area->h; y++) {
			for (x = 0; x < area->w; x++) {
				buffer[y*area->w + x] = rgb_buf[y*area->w*3 + x*3 + color_offset]; 
			}
		}
		free(rgb_buf);
		break;
	default:
		siril_log_message("This type of Bayer pattern is not handled yet.\n");
		return -1;
	}

	return 0;
}

int ser_write_frame_from_fit(struct ser_struct *ser_file, fits *fit, int frame_no) {
	int frame_size, pixel, plane, dest;
	int ret, retval = 0;
	off_t offset;
	BYTE *data8 = NULL;			// for 8-bit files
	WORD *data16 = NULL;		// for 16-bit files

	if (!ser_file || ser_file->fd <= 0 || !fit)
		return -1;
	if (ser_file->number_of_planes == 0) {
		// adding first frame of a new sequence, use it to populate the header
		ser_header_from_fit(ser_file, fit);
	}
	if (fit->rx != ser_file->image_width || fit->ry != ser_file->image_height) {
		siril_log_message("Trying to add an image of different size in a SER\n");
		return 1;
	}

	fits_flip_top_to_bottom(fit);
	frame_size = ser_file->image_width * ser_file->image_height *
		ser_file->number_of_planes;

	offset = SER_HEADER_LEN	+ (off_t) frame_size *
			(off_t) ser_file->byte_pixel_depth * (off_t) frame_no;

	if (ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_8)
		data8 = malloc(frame_size * ser_file->byte_pixel_depth);
	else
		data16 = malloc(frame_size * ser_file->byte_pixel_depth);

	for (plane = 0; plane < ser_file->number_of_planes; plane++) {
		dest = plane;
		for (pixel = 0; pixel < ser_file->image_width * ser_file->image_height;
				pixel++) {
			if (ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_8)
				data8[dest] = (BYTE)(fit->pdata[plane][pixel]);
			else data16[dest] = (fit->pdata[plane][pixel] >> 8 | fit->pdata[plane][pixel] << 8);
			dest += ser_file->number_of_planes;
		}
	}

#ifdef _OPENMP
	omp_set_lock(&ser_file->fd_lock);
#endif
	if ((off_t)-1 == lseek(ser_file->fd, offset, SEEK_SET)) {
#ifdef _OPENMP
		omp_unset_lock(&ser_file->fd_lock);
#endif
		perror("seek");
		retval = -1;
		goto free_and_quit;
	}

	if (ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_8) {
		ret = write(ser_file->fd, data8, frame_size * ser_file->byte_pixel_depth);
#ifdef _OPENMP
		omp_unset_lock(&ser_file->fd_lock);
#endif
		if (ret != frame_size * ser_file->byte_pixel_depth) {
			perror("write image in SER");
			retval = 1;
			goto free_and_quit;
		}
	} else {
		ret = write(ser_file->fd, data16, frame_size * ser_file->byte_pixel_depth);
#ifdef _OPENMP
		omp_unset_lock(&ser_file->fd_lock);
#endif
		if (ret != frame_size * ser_file->byte_pixel_depth) {
			perror("write image in SER");
			retval = 1;
			goto free_and_quit;
		}
	}

#ifdef _OPENMP
#pragma omp atomic
#endif
	ser_file->frame_count++;
	retval = 0;

free_and_quit:
	if (data8) free(data8);
	if (data16) free(data16);
	return retval;
}

void set_combo_box_bayer_pattern(ser_color pattern) {
	GtkComboBox *combo = GTK_COMBO_BOX(lookup_widget("comboBayer_pattern"));
	int entry;

	switch (pattern) {
	default:
	case SER_BAYER_RGGB:
		entry = 0;
		break;
	case SER_BAYER_BGGR:
		entry = 1;
		break;
	case SER_BAYER_GBRG:
		entry = 2;
		break;
	case SER_BAYER_GRBG:
		entry = 3;
		break;
	}
	gtk_combo_box_set_active(combo, entry);
}

/* once a buffer (data) has been acquired from the file, with frame_size pixels
 * read in it, depending on ser_file's endianess and pixel depth, data is
 * reorganized to match Siril's data format. */
void ser_manage_endianess_and_depth(struct ser_struct *ser_file, WORD *data, int frame_size) {
	WORD pixel;
	int i;
	if (ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_8) {
		// inline conversion to 16 bit
		for (i = frame_size - 1; i >= 0; i--)
			data[i] = (WORD) (((BYTE*)data)[i]);
	} else if (ser_file->little_endian) {	// TODO check if it is needed for big endian
		// inline conversion to big endian
		for (i = frame_size - 1; i >= 0; i--) {
			pixel = data[i];
			pixel = (pixel >> 8) | (pixel << 8);
			data[i] = pixel;
		}
	}
}

