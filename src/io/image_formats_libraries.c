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
*/

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#ifdef HAVE_LIBTIFF
#define uint64 uint64_hack_
#define int64 int64_hack_
#include <tiffio.h>
#undef uint64
#undef int64
#endif
#ifdef HAVE_LIBJPEG
#include <jpeglib.h>
#endif
#ifdef HAVE_LIBPNG
#include <png.h>
#include <setjmp.h>
#endif
#ifdef HAVE_LIBRAW
#include <libraw/libraw.h>
#endif
#ifdef HAVE_LIBGIF
#include <gif_lib.h>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "gui/callbacks.h"
#include "single_image.h"

/********************* TIFF IMPORT AND EXPORT *********************/

#ifdef HAVE_LIBTIFF

/* reads a TIFF file and stores it in the fits argument.
 * If file loading fails, the argument is untouched.
 */
int readtif(const char *name, fits *fit) {
	char *msg;
	int retval = 0;
	uint32 height, width, npixels;
	uint16 nbits, nsamples, color;
	WORD *data;
	
	TIFF* tif = TIFFOpen(name, "r");
	if (!tif) {
		siril_log_message(_("Could not open the TIFF file %s\n"), name);
		return -1;
	}
	data = NULL;
	
	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
	TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &nsamples);	
	TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &nbits);
	TIFFGetField(tif, TIFFTAG_PHOTOMETRIC, &color);
	TIFFGetField(tif, TIFFTAG_MINSAMPLEVALUE, &(fit->lo));
	TIFFGetField(tif, TIFFTAG_MAXSAMPLEVALUE, &(fit->hi));

	npixels = width * height;

	switch(nbits){
		case 8:
			/* High level functions in readtif8bits: should read every 8-bit TIFF file */
			retval = readtif8bits(tif, width, height, nsamples, &data);
			break;

		case 16:
			/* Lower level functions in readtifstrip: 
			 * maybe some particular TIFF files could not be read.
			 * No such files are known yet */
			retval = readtifstrip(tif, width, height, nsamples, &data);
			break;

		default :
			msg = siril_log_message(_("Siril only works with 8/16-bit TIFF format.\n"));
			show_dialog(msg, _("Warning"), "gtk-dialog-warning");
			retval = -1;
	}
	TIFFClose(tif);
	if (retval <= 0) {
		if (data)
			free(data);
	} else if (data != NULL) {
		clearfits(fit);
		fit->rx = width;
		fit->ry = height;
		fit->naxes[0] = width;
		fit->naxes[1] = height;
		fit->data = data;
		fit->binning_x=fit->binning_y=1;
		if (nsamples == 1) {
			fit->naxes[2] = nsamples;
			fit->naxis = 2;
			fit->pdata[RLAYER]=fit->data;
			fit->pdata[GLAYER]=fit->data;
			fit->pdata[BLAYER]=fit->data;
		} else {
			fit->naxes[2] = 3;
			fit->naxis = 3;
			fit->pdata[RLAYER]=fit->data;
			fit->pdata[GLAYER]=fit->data + npixels;
			fit->pdata[BLAYER]=fit->data + npixels * 2;
		}
		fit->bitpix = (nbits == 8) ? BYTE_IMG : USHORT_IMG;
		retval = nsamples;
	}
	if (nbits==16) mirrorx(fit, FALSE);
	siril_log_message(_("Reading TIFF: %d-bit file %s, %ld layer(s), %ux%u pixels\n"),
						nbits, name, fit->naxes[2], fit->rx, fit->ry);
	return retval;
}

int readtifstrip(TIFF* tif, uint32 width, uint32 height, uint16 nsamples, WORD **data) {
	unsigned int npixels;
	int  i, j, scanline, retval=nsamples;
	WORD *buf;
	uint32 rowsperstrip;
	uint16 config;
	u_long nrow, row;
	char *msg;

	TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
	TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip);
	   
	npixels = width * height;
	*data = malloc(npixels * sizeof(WORD) * nsamples);
	WORD *gbuf[3] = {*data, *data, *data};	
	if (nsamples == 4) {
		msg = siril_log_message(_("Alpha channel is ignored.\n"));
		show_dialog(msg, _("Warning"), "gtk-dialog-warning");
	}
	if ((nsamples == 3) || (nsamples == 4)) {
		gbuf[1] = *data + npixels;
		gbuf[2] = *data + npixels * 2;
	}

	scanline = TIFFScanlineSize(tif);
	buf = _TIFFmalloc(TIFFStripSize(tif));
	if (!buf) return -1;
	for (row = 0; row < height; row += rowsperstrip){
		nrow = (row + rowsperstrip > height ? height - row : rowsperstrip);
		switch(config){
			case PLANARCONFIG_CONTIG:
				if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, row, 0), buf, nrow*scanline) < 0){
					msg = siril_log_message(_("An unexpected error was encountered while trying to read the file.\n"));
					show_dialog(msg, _("Error"), "gtk-dialog-error");
					retval = -1;
					break;
				}
				for (i = 0; i < width*nrow; i++) {
					*gbuf[RLAYER]++ = buf[i*nsamples+0];
					if ((nsamples == 3) || (nsamples == 4)) {
						*gbuf[GLAYER]++ = buf[i*nsamples+1];
						*gbuf[BLAYER]++ = buf[i*nsamples+2];
					}
				}
				break;
			case PLANARCONFIG_SEPARATE:	
				if (nsamples >= 3)		//don't need to read the alpha
					nsamples = 3;
				for (j=0; j<nsamples; j++){	//loop on the layer
					if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, row, j), buf, nrow*scanline) < 0){
						msg = siril_log_message(_("An unexpected error was encountered while trying to read the file.\n"));
						show_dialog(msg, _("Error"), "gtk-dialog-error");
						retval = -1;
						break;
					}
					for (i = 0; i < width*nrow; i++)
						*gbuf[j]++ = buf[i];
				}
				break;
			default:
				msg = siril_log_message(_("Unknown TIFF file.\n"));
				show_dialog(msg, _("Error"), "gtk-dialog-error");
				retval = -1;
		}
	}
	_TIFFfree(buf);
	return retval;
}
int readtif8bits(TIFF* tif, uint32 width, uint32 height, uint16 nsamples, WORD **data) {
	uint32 npixels;
	int retval = nsamples;
	char *msg;
    
	npixels = width * height;
	*data = malloc(npixels * sizeof(WORD) * nsamples);
	WORD *gbuf[3] = {*data, *data, *data};
	if (nsamples == 4) {
		siril_log_message(_("Alpha channel is ignored.\n"));
	}
	if ((nsamples == 3) || (nsamples == 4)) {
		gbuf[1] = *data + npixels;
		gbuf[2] = *data + npixels * 2;
	}

	/* get the data */
	uint32* raster = (uint32*)_TIFFmalloc(npixels*sizeof(uint32));
	if (raster != NULL) {
		if (TIFFReadRGBAImage(tif, width, height, raster, 0)) {
			int i, j;
			for (j = 0; j < height; j++) {
				int istart = j * width;
				for (i = 0; i < width; i++) {
					*gbuf[RLAYER]++ = (WORD)TIFFGetR(raster[istart + i]);
					if ((nsamples == 3) || (nsamples == 4)) {
						*gbuf[GLAYER]++ = (WORD)TIFFGetG(raster[istart + i]);
						*gbuf[BLAYER]++ = (WORD)TIFFGetB(raster[istart + i]);
					}
				}
			}
		}
		else {
			msg = siril_log_message(_("An unexpected error was encountered while trying to read the file.\n"));
			show_dialog(msg, _("Error"), "gtk-dialog-error");
			retval = -1;
		}
		_TIFFfree(raster);
	}
	else retval = -1;
	return retval;
}

/*** This function save the current image into a uncompressed 8- or 16-bit file *************/

int savetif(const char *name, fits *fit, uint16 bitspersample){
	int retval = 0;
	char *msg;
	unsigned char *buf8;
	WORD *buf16;
	uint32 width, height, row, col, n;
	uint16 nsamples;

	mirrorx(fit, FALSE);

	TIFF* tif = TIFFOpen(name, "w");

	if (tif == NULL) {
		msg = siril_log_message(_("Siril cannot create TIFF file.\n"));
		show_dialog(msg, _("Error"), "gtk-dialog-error");
		return 1;
	}
	nsamples	=(uint16)fit->naxes[2];
	width		=(uint32)fit->rx;
	height		=(uint32)fit->ry;	
	
	/*******************************************************************
	 * If the user saves a tif from the graphical menu, he can set
	 * the Description and the Copyright of the Image 
	 ******************************************************************/
	gchar *img_desc=NULL, *img_copy=NULL;
	GtkTextView *description = GTK_TEXT_VIEW(lookup_widget("Description_txt"));
	GtkTextView *copyright = GTK_TEXT_VIEW(lookup_widget("Copyright_txt"));
	GtkTextBuffer *desbuf = gtk_text_view_get_buffer(description);
	GtkTextBuffer *copybuf = gtk_text_view_get_buffer(copyright);
	GtkTextIter itDebut;
	GtkTextIter itFin;

	gtk_text_buffer_get_start_iter(desbuf, &itDebut);
	gtk_text_buffer_get_end_iter(desbuf, &itFin);
	if (desbuf)
		img_desc=gtk_text_buffer_get_text(desbuf, &itDebut, &itFin, TRUE);
	gtk_text_buffer_get_bounds (desbuf, &itDebut, &itFin);
	gtk_text_buffer_delete (desbuf, &itDebut, &itFin);
	gtk_text_buffer_get_start_iter(copybuf, &itDebut);
	gtk_text_buffer_get_end_iter(copybuf, &itFin);
	if (copybuf)
		img_copy=gtk_text_buffer_get_text(copybuf, &itDebut, &itFin, TRUE);
	gtk_text_buffer_get_bounds (copybuf, &itDebut, &itFin);
	gtk_text_buffer_delete (copybuf, &itDebut, &itFin);
		
	/*******************************************************************/

	/* TIFF TAG FIELD */
	TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bitspersample);
	TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
	TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
	TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tif, -1));
	TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, nsamples);
	TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(tif, TIFFTAG_IMAGEDESCRIPTION, img_desc);
	TIFFSetField(tif, TIFFTAG_COPYRIGHT, img_copy);
	TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(tif, TIFFTAG_MINSAMPLEVALUE, fit->mini);
	TIFFSetField(tif, TIFFTAG_MAXSAMPLEVALUE, fit->maxi);
	TIFFSetField(tif, TIFFTAG_SOFTWARE, PACKAGE " v" VERSION);
	
	if (nsamples == 1)
		TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	else if (nsamples == 3)
		TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
	else {
		TIFFClose(tif);
		msg = siril_log_message(_("TIFF file has unexpected number of channels (not 1 or 3).\n"));
		show_dialog(msg, _("Error"), "gtk-dialog-error");
		return 1;
	}

	WORD *gbuf[3]={fit->pdata[RLAYER],fit->pdata[GLAYER],fit->pdata[BLAYER]};  
	
	switch (bitspersample){
		case 8:
			buf8 = _TIFFmalloc(width * sizeof(unsigned char) * nsamples);
			for (row = 0; row < height; row++){
				for (col = 0; col < width; col++){
					for (n=0; n < nsamples; n++) {
						/* UCHAR_MAX / USHRT_MAX is constant, it's 1/255
						 * This operation should be speed-up by doing a shift */
						buf8[col*nsamples+n]= (*gbuf[n]++) * UCHAR_MAX_DOUBLE/USHRT_MAX_DOUBLE;
						//~ buf8[col*nsamples+n] = *gbuf[n]++ >> 8; // bad shift
					}
				}
				TIFFWriteScanline(tif, buf8, row, 0);
			}
		_TIFFfree(buf8);
			break;
		case 16: 
			buf16 = _TIFFmalloc(width * sizeof(WORD) * nsamples);
			for (row = 0; row < height; row++){
				for (col = 0; col < width; col++){
					for (n=0; n < nsamples; n++)
						buf16[col*nsamples+n]=(*gbuf[n]++);
				}
				TIFFWriteScanline(tif, buf16, row, 0);
			}
		_TIFFfree(buf16);
		break;
		default:		// Should not happen
			retval = 1;
	}
	TIFFClose(tif);
	mirrorx(fit, FALSE);
	siril_log_message(_("Saving TIFF: %d-bit file %s, %ld layer(s), %ux%u pixels\n"),
						bitspersample, name, fit->naxes[2], fit->rx, fit->ry);
	return retval;
}
#endif	// HAVE_LIBTIFF


/********************* JPEG IMPORT AND EXPORT *********************/

#ifdef HAVE_LIBJPEG
int readjpg(const char* name, fits *fit){
	int i, npixels;
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	WORD *data;
	FILE *f;
	JSAMPARRAY pJpegBuffer;
	int row_stride;
	if ((f = fopen(name, "rb")) == NULL){
		char *msg = siril_log_message(_("Sorry but Siril cannot open the file: %s.\n"), name);
		show_dialog(msg, _("Error"), "gtk-dialog-error");
		return -1;
	}
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, f);
	jpeg_read_header(&cinfo, TRUE);
	jpeg_start_decompress(&cinfo);
	npixels=cinfo.output_width*cinfo.output_height;
	data = malloc(npixels * sizeof(WORD) * 3);
	WORD *buf[3] = {data, data + npixels, data + npixels * 2};
	row_stride = cinfo.output_width * cinfo.output_components ;
	pJpegBuffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
	while (cinfo.output_scanline < cinfo.output_height) {
		jpeg_read_scanlines(&cinfo, pJpegBuffer, 1);
		for (i=0;i<cinfo.output_width;i++) {
			*buf[RLAYER]++ = pJpegBuffer[0][cinfo.output_components*i+0];
			*buf[GLAYER]++ = pJpegBuffer[0][cinfo.output_components*i+1];
			*buf[BLAYER]++ = pJpegBuffer[0][cinfo.output_components*i+2];
		}
	}
	fclose(f);
	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	if (data != NULL){
		clearfits(fit);	
		fit->bitpix = BYTE_IMG;
		if (cinfo.output_components==1)
			fit->naxis = 2;
		else
			fit->naxis = 3;
		fit->rx = cinfo.output_width;
		fit->ry = cinfo.output_height;
		fit->naxes[0] = cinfo.output_width;
		fit->naxes[1] = cinfo.output_height;
		fit->naxes[2] = cinfo.output_components;
		fit->data = data;
		fit->pdata[RLAYER]=fit->data;
		fit->pdata[GLAYER]=fit->data + npixels;
		fit->pdata[BLAYER]=fit->data + npixels * 2;
		fit->binning_x=fit->binning_y=1;
	}
	mirrorx(fit, FALSE);
	siril_log_message(_("Reading JPG: file %s, %ld layer(s), %ux%u pixels\n"),
						name, fit->naxes[2], fit->rx, fit->ry);
	return cinfo.output_components;
}

int savejpg(char *name, fits *fit, int quality){
	FILE *f;
	int i, j;
	unsigned char red, blue, green;
	struct jpeg_compress_struct cinfo;    // Basic info for JPEG properties.
	struct jpeg_error_mgr jerr;           // In case of error.
	JSAMPROW row_pointer[1];              // Pointer to JSAMPLE row[s].
	int row_stride;                       // Physical row width in image buffer.

	//## ALLOCATE AND INITIALIZE JPEG COMPRESSION OBJECT
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);

	//## OPEN FILE FOR DATA DESTINATION:
	if ((f = fopen(name, "wb")) == NULL) {
		char *msg = siril_log_message(_("Siril cannot create JPG file.\n"));
		show_dialog(msg, _("Error"), "gtk-dialog-error");
		return 1;
	}
	jpeg_stdio_dest(&cinfo, f);

	//## SET PARAMETERS FOR COMPRESSION:
	cinfo.image_width  = fit->rx;   // |-- Image width and height in pixels.
	cinfo.image_height = fit->ry;   // |
	cinfo.input_components = 3;     // Number of color components per pixel.
	cinfo.in_color_space = JCS_RGB; // Colorspace of input image as RGB.

	unsigned char *gbuf[3]={
		com.graybuf[RLAYER]+cinfo.image_width*(cinfo.image_height-1)*4,
		com.graybuf[GLAYER]+cinfo.image_width*(cinfo.image_height-1)*4,
		com.graybuf[BLAYER]+cinfo.image_width*(cinfo.image_height-1)*4 };
	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, quality, TRUE);

	//## CREATE IMAGE BUFFER TO WRITE FROM AND MODIFY THE IMAGE TO LOOK LIKE CHECKERBOARD:
	unsigned char *image_buffer = NULL;
	image_buffer = (unsigned char*)malloc(cinfo.image_width*cinfo.image_height*cinfo.num_components);
	for (i=(cinfo.image_height-1);i>=0;i--){
		for (j=0;j<cinfo.image_width;j++){
			int pixelIdx = ((i*cinfo.image_width)+j) * cinfo.input_components;
			red = *gbuf[RLAYER];
			gbuf[RLAYER] += 4;
			if (fit->naxes[2]==3){
				green   = *gbuf[GLAYER];
				blue    = *gbuf[BLAYER];
				gbuf[GLAYER] += 4;
				gbuf[BLAYER] += 4;
			} else {
				green   = red;
				blue    = red;
			}
			image_buffer[pixelIdx+0] = red;         // r |-- Set r,g,b components to
			image_buffer[pixelIdx+1] = green;       // g |   make this pixel 
			image_buffer[pixelIdx+2] = blue;        // b |
		}
		if (fit->naxes[2]==3){
			gbuf[GLAYER] -= cinfo.image_width * 8;
			gbuf[BLAYER] -= cinfo.image_width * 8;
		}
		gbuf[RLAYER] -= cinfo.image_width * 8;
	}
	//## START COMPRESSION:
	jpeg_start_compress(&cinfo, TRUE);
	row_stride = cinfo.image_width * 3;        // JSAMPLEs per row in image_buffer

	while (cinfo.next_scanline < cinfo.image_height)
	{
		row_pointer[0] = &image_buffer[cinfo.next_scanline * row_stride];
		(void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}
	// NOTE: jpeg_write_scanlines expects an array of pointers to scanlines.
	//       Here the array is only one element long, but you could pass
	//       more than one scanline at a time if that's more convenient.

	//## FINISH COMPRESSION AND CLOSE FILE:
	jpeg_finish_compress(&cinfo);
	fclose(f);
	jpeg_destroy_compress(&cinfo);
	free(image_buffer);
	siril_log_message(_("Saving JPG: file %s, quality=%d%%, %ld layer(s), %ux%u pixels\n"),
						name, quality, fit->naxes[2], fit->rx, fit->ry);
	return 0;
}
#endif	// HAVE_LIBJPEG


/********************* PNG IMPORT *********************/

#ifdef HAVE_LIBPNG
/* reads a PNG file and stores it in the fits argument.
 */
int readpng(const char *name, fits* fit) {
	int width, height, x, y, npixels, nbplanes;
	png_byte color_type;
	png_byte bit_depth;
	png_bytep *row_pointers;
	WORD *data;
	FILE *f;

	png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if(!png) {
		char *msg = siril_log_message(_("Sorry but Siril cannot open the file: %s.\n"), name);
		show_dialog(msg, _("Error"), "gtk-dialog-error");
		return -1;
	}

	png_infop info = png_create_info_struct(png);
	if(!info) return -1;

	if(setjmp(png_jmpbuf(png))) return -1;

	f = fopen(name, "rb");
	png_init_io(png, f);

	png_read_info(png, info);

	width = png_get_image_width(png, info);
	height = png_get_image_height(png, info);
	npixels = width * height;
	color_type = png_get_color_type(png, info);
	bit_depth = png_get_bit_depth(png, info);

	data = malloc(npixels * sizeof(WORD) * 3);
	WORD *buf[3] = {data, data + npixels, data + npixels * 2};

	if(color_type == PNG_COLOR_TYPE_PALETTE)
		png_set_palette_to_rgb(png);

	// PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth.
	if(color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
		png_set_expand_gray_1_2_4_to_8(png);

	if(png_get_valid(png, info, PNG_INFO_tRNS))
		png_set_tRNS_to_alpha(png);

	// These color_type don't have an alpha channel then fill it with 0xff.
	if(color_type == PNG_COLOR_TYPE_RGB || color_type == PNG_COLOR_TYPE_GRAY ||
			color_type == PNG_COLOR_TYPE_PALETTE)
		png_set_filler(png, 0xFF, PNG_FILLER_AFTER);

	if(color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
		png_set_gray_to_rgb(png);

	png_read_update_info(png, info);

	row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * height);
	for(y = 0; y < height; y++) {
		row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png,info));
	}

	png_read_image(png, row_pointers);

	fclose(f);
	
	if (bit_depth == 16){			//in 16-bit: it is stored as RRGGBB
		for (y=0; y<height; y++) {
			png_byte* row = row_pointers[y];
			for (x=0; x<width; x++) {
				png_byte* ptr = &(row[x*8]);
				*buf[RLAYER]++ = ptr[0]*(UCHAR_MAX+1)+ptr[1];
				*buf[GLAYER]++ = ptr[2]*(UCHAR_MAX+1)+ptr[3];
				*buf[BLAYER]++ = ptr[4]*(UCHAR_MAX+1)+ptr[5];	
			}
		}
	}
	else{
		for (y=0; y<height; y++) {
			png_byte* row = row_pointers[y];
			for (x=0; x<width; x++) {		
				png_byte* ptr = &(row[x*4]);
				*buf[RLAYER]++ = ptr[0];
				*buf[GLAYER]++ = ptr[1];
				*buf[BLAYER]++ = ptr[2];
			}
		}
		fit->bitpix = BYTE_IMG;
	}	
	// We define the number of channel we have
	switch(color_type){
	case PNG_COLOR_TYPE_RGB_ALPHA:
	case PNG_COLOR_TYPE_RGB:
		nbplanes=3;
		break;
	case PNG_COLOR_TYPE_GRAY_ALPHA:
	case PNG_COLOR_TYPE_GRAY:
		nbplanes=1;
		break;
	default:
		nbplanes=0;	
	}
	
	// free allocated memory
	for (y=0; y<height; y++)
		free(row_pointers[y]);
	free(row_pointers);
	
	if (data != NULL){
		clearfits(fit);
		fit->rx = width;
		fit->ry = height;
		fit->naxes[0] = width;
		fit->naxes[1] = height;
		fit->naxes[2] = nbplanes;
		if (nbplanes==1)
			fit->naxis = 2;
		else
			fit->naxis = 3;
		fit->bitpix = (bit_depth == 8) ? BYTE_IMG : USHORT_IMG;
		fit->data = data;
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->data + npixels;
		fit->pdata[BLAYER] = fit->data + npixels * 2;
		fit->binning_x = fit->binning_y = 1;

	}
	mirrorx(fit, FALSE);
	siril_log_message(_("Reading PNG: %d-bit file %s, %ld layer(s), %ux%u pixels\n"),
						bit_depth, name, fit->naxes[2], fit->rx, fit->ry);

	return nbplanes;
}
#endif	// HAVE_LIBPNG

/********************* RAW IMPORT *********************/
#ifdef HAVE_LIBRAW

int readraw(const char *name, fits *fit) {
	ushort width, height;
	int npixels, i;
	WORD *data=NULL;
	libraw_data_t *raw = libraw_init(0);
	libraw_processed_image_t *image = NULL;
	int ret = libraw_open_file(raw, name);
	if (ret) {
		siril_log_message(_("Error in libraw %s\n"), libraw_strerror(ret));
		libraw_recycle(raw);
		libraw_close(raw);	
		return -1;
	}
	
	if (raw->other.shutter > 1.0)
		siril_log_message(_("Decoding %s %s file (ISO=%g, Exposure=%gs)\n"),
				raw->idata.make, raw->idata.model, raw->other.iso_speed, raw->other.shutter);
	else
		siril_log_message(_("Decoding %s %s file (ISO=%g, Exposure=1/%gs)\n"),
				raw->idata.make, raw->idata.model, raw->other.iso_speed, 1/raw->other.shutter);
		
	raw->params.output_bps = 16;						/* 16-bits files                           */
	raw->params.four_color_rgb = 0;						/* If == 1, interpolate RGB as four colors.*/
	raw->params.no_auto_bright = 1;						/* no auto_bright                          */
	raw->params.gamm[0] = 1./com.raw_set.gamm[0];		/* Gamma curve set by the user             */
	raw->params.gamm[1] = com.raw_set.gamm[1];	                                         
	raw->params.bright = com.raw_set.bright;			/* Brightness                              */
	raw->params.user_flip = 0;							/* no flip                                 */
	raw->params.use_camera_wb = com.raw_set.use_camera_wb;
	raw->params.use_auto_wb = com.raw_set.use_auto_wb;
	if (com.raw_set.user_black==1)
		raw->params.user_black = 0;						/* user black level equivalent to dcraw -k 0 */
	raw->params.output_color = 0;						/* output colorspace, 0=raw, 1=sRGB, 2=Adobe, 3=Wide, 4=ProPhoto, 5=XYZ*/
		
	if (!(com.raw_set.auto_mul)) {		/* 4 multipliers (r,g,b,g) of the user's white balance.    */
		raw->params.user_mul[0] = (float)com.raw_set.mul[0];
		raw->params.user_mul[1] = raw->params.user_mul[3] = 1.0f;
		raw->params.user_mul[2] = (float)com.raw_set.mul[2];
		siril_log_message(_("Daylight multipliers: %f, %f, %f\n"),
				raw->params.user_mul[0], raw->params.user_mul[1], raw->params.user_mul[2]);
	}
	else {
		float mul[4];					/* 3 multipliers (r,g,b) from the camera white balance.  */
		mul[0] = raw->color.pre_mul[0]/raw->color.pre_mul[1];
		mul[1] = 1.0; /* raw->color.pre_mul[1]/raw->color.pre_mul[1]; */
		mul[2] = raw->color.pre_mul[2]/raw->color.pre_mul[1];
		mul[3] = raw->color.pre_mul[3]/raw->color.pre_mul[1];
		siril_log_message(_("Daylight multipliers: %f, %f, %f\n"), mul[0], mul[1], mul[2]);
	}

	switch(com.raw_set.user_qual) {		/* Set interpolation                                        */
		case 0:		/* bilinear interpolaton */
			raw->params.user_qual = 0;
			siril_log_message(_("Bilinear interpolation...\n"));
			break;
		case 2:		/* VNG interpolaton */
			raw->params.user_qual = 1;
			siril_log_message(_("VNG interpolation...\n"));
			break;
		case 3:		/* PPG interpolaton */
			raw->params.user_qual = 2;
			siril_log_message(_("PPG interpolation...\n"));
			break;
		default:
		case 1:		/* AHD interpolaton */
			raw->params.user_qual = 3;
			siril_log_message(_("AHD interpolation...\n"));
			break;
	}
	
	width = raw->sizes.iwidth;
	height = raw->sizes.iheight;
	npixels = width * height;
	
	data = malloc(npixels * sizeof(WORD) * 3);
	if (!data) return -1;
	WORD *buf[3] = {data, data + npixels, data + npixels * 2};
	ret = libraw_unpack(raw);
	if (ret) {
		printf("Error in libraw %s\n", libraw_strerror(ret));
		free(data);
		libraw_recycle(raw);
		libraw_close(raw);
		return -1;
	}
	
	ret = libraw_dcraw_process(raw);
	if (ret) {
		printf("Error in libraw %s\n", libraw_strerror(ret));
		free(data);
		libraw_recycle(raw);
		libraw_close(raw);
		return -1;
	}
	image = libraw_dcraw_make_mem_image(raw, &ret);
	if (ret) {
		printf("Error in libraw %s\n", libraw_strerror(ret));
		free(data);
		libraw_dcraw_clear_mem(image);
		libraw_recycle(raw);
		libraw_close(raw);	
		return -1;
	}

	int nbplanes = image->colors;
	if (nbplanes!=3) {
		free(data);
		libraw_dcraw_clear_mem(image);
		libraw_recycle(raw);
		libraw_close(raw);	
		return -1;
	}
	// only for 16-bits because of endianness. Are there 8-bits RAW ???

	for (i=0; i<image->data_size; i+=6) {
		*buf[RLAYER]++ = (image->data[i+0]) + (image->data[i+1] << 8);
		*buf[GLAYER]++ = (image->data[i+2]) + (image->data[i+3] << 8);
		*buf[BLAYER]++ = (image->data[i+4]) + (image->data[i+5] << 8);
	}

	libraw_dcraw_clear_mem(image);
	libraw_recycle(raw);
	libraw_close(raw);	
	
	if (data != NULL) {
		clearfits(fit);
		fit->bitpix = USHORT_IMG;
		fit->rx = width;
		fit->ry = height;
		fit->naxes[0] = width;
		fit->naxes[1] = height;
		fit->naxes[2] = nbplanes;
		if (nbplanes==1)
			fit->naxis = 2;
		else
			fit->naxis = 3;	
		fit->data = data;
		fit->pdata[RLAYER]=fit->data;
		fit->pdata[GLAYER]=fit->data + npixels;
		fit->pdata[BLAYER]=fit->data + npixels * 2;
		fit->binning_x=fit->binning_y=1;
		if (raw->other.focal_len>0.) fit->focal_length=raw->other.focal_len;
		if (raw->other.iso_speed>0.) fit->iso_speed=raw->other.iso_speed;
		if (raw->other.shutter>0.) fit->exposure=raw->other.shutter;
		if (raw->other.aperture>0.) fit->aperture=raw->other.aperture;
		snprintf(fit->instrume, FLEN_VALUE, "%s %s", raw->idata.make, raw->idata.model);
	}
	return nbplanes;
}

int readraw_in_cfa(const char *name, fits *fit) {
	libraw_data_t *raw = libraw_init(0);
	unsigned int i, col, row;
	ushort raw_width, width, height, left_margin, top_margin;
	int npixels;
	WORD *data=NULL;
	int ret = libraw_open_file(raw, name);
	
	if (ret) {
		printf("Error in libraw %s\n", libraw_strerror(ret));
		return -1;
	}
	
	ret = libraw_unpack(raw);
	if (ret) {
		printf("Error in libraw %s\n", libraw_strerror(ret));
		return -1;
	}
	
	/* This test checks if raw data exist. Sometimes it doesn't. This is
	 * the case for DNG built from lightroom for example */
	if ((void*)raw->rawdata.raw_image != raw->rawdata.raw_alloc || raw->rawdata.raw_alloc == 0x00){
		siril_log_message(_("Siril cannot open this file in CFA mode (no data available). Try to switch into RGB.\n"));
		return -1;
	}
		
	raw->params.user_flip = 0;				/* no flip                                 */
	raw->params.output_color = 0;			/* output colorspace, 0=raw, 1=sRGB, 2=Adobe, 3=Wide, 4=ProPhoto, 5=XYZ*/

	raw_width = raw->sizes.raw_width;
	width = raw->sizes.width;
	height = raw->sizes.height;
	npixels = width * height;	
	left_margin = raw->rawdata.sizes.left_margin;
	top_margin = raw->rawdata.sizes.top_margin;
	
	if (raw->other.shutter > 1.0)
		siril_log_message(_("Decoding %s %s file (ISO=%g, Exposure=%gs)\n"),
				raw->idata.make, raw->idata.model, raw->other.iso_speed, raw->other.shutter);
	else
		siril_log_message(_("Decoding %s %s file (ISO=%g, Exposure=1/%gs)\n"),
				raw->idata.make, raw->idata.model, raw->other.iso_speed, 1/raw->other.shutter);

	data = (WORD*) calloc(1, npixels * sizeof(WORD));
	if (!data) {
		return -1;
	}
	WORD *buf = data;
	
	int offset = raw_width * top_margin + left_margin;	
	i = 0;
	for (row = 0; row < height; row++){
		for (col = 0; col < width; col++){
			buf[i++] = (WORD)raw->rawdata.raw_image[offset + col + (raw_width * row)] ;
		}
	}

	if (data != NULL) {
		clearfits(fit);
		fit->bitpix = USHORT_IMG;
		fit->rx = (unsigned int) width;
		fit->ry = (unsigned int) height;
		fit->naxes[0] = (long) width;
		fit->naxes[1] = (long) height;
		fit->naxes[2] = 1;
		fit->naxis = 2;
		fit->data = data;
		fit->pdata[RLAYER]=fit->data;
		fit->pdata[GLAYER]=fit->data;
		fit->pdata[BLAYER]=fit->data;
		fit->binning_x=fit->binning_y=1;
		if (raw->other.focal_len>0.) fit->focal_length=raw->other.focal_len;
		if (raw->other.iso_speed>0.) fit->iso_speed=raw->other.iso_speed;
		if (raw->other.shutter>0.) fit->exposure=raw->other.shutter;
		if (raw->other.aperture>0.) fit->aperture=raw->other.aperture;
		snprintf(fit->instrume, FLEN_VALUE, "%s %s", raw->idata.make, raw->idata.model);
	}

	libraw_recycle(raw);
	libraw_close(raw);		
	return 1;
}

int open_raw_files(const char *name, fits *fit, int type) {
	int retvalue = 1;
	
	switch(type){
		default:
		case 0:
			retvalue = readraw(name, fit);
			break;
		case 1:
			retvalue = readraw_in_cfa(name, fit);
			break;
	}
	mirrorx(fit, FALSE);
	siril_log_message(_("Reading RAW: file %s, %ld layer(s), %ux%u pixels\n"),
						name, fit->naxes[2], fit->rx, fit->ry);		
	return retvalue;
}
#endif


/********************* GIF EXPORT *********************/
#ifdef HAVE_LIBGIF

/* save a fit as gif.
 * If anim is true, if gif is nul, it is set to the created GIF handle. If it's
 * not nul it's used as is and a new image is added to the GIF.
 * If an error occurs, everything is freed and the file may be closed, in that
 * case, gif is set to NULL.
 * No check is done on filename actually corresponding to gif.
 * if anim, gif must be finalized with gifclose(gif).
 * if positive, delay and loop_count are used for animation. Delay is in centiseconds.
 *
 * palette code inspired by https://github.com/schani/animatedgif/blob/master/writegif.c
 */
int savegif(const char *filename, fits *fit, int anim, GifFileType **gif, int delay, int loop_count) {
	/* convert the fits to RGB 8-bit buffers, using the min/max scale */
	BYTE *rgb[3];
	BYTE map[USHRT_MAX + 1];
	GifByteType *gif_buffer;
	WORD tmp_pixel_value, hi, lo;
	int i, nb_pixels;
	float pente;

	if (anim && gif == NULL) {
		printf("Gif should not be NULL in anim mode\n");
		return 1;
	}

	nb_pixels = fit->rx * fit->ry;
	if (fit->naxes[2] == 1) {
		/*  gif_buffer contains values of the palette, which in our case of
		 *  linear grey palette are the same as the pixel values */
		gif_buffer = malloc(nb_pixels);
		if (!gif_buffer) {
			printf("Failed to allocate memory required, aborted.\n");
			return 1;
		}
	} else {
		assert(fit->naxes[2] == 3);
		gif_buffer = NULL;

		rgb[0] = malloc(nb_pixels);
		if (rgb[0] == NULL) {
			printf("Failed to allocate memory required, aborted.\n");
			if (anim) *gif = NULL;
			return 1;
		}
		rgb[1] = malloc(nb_pixels);
		if (rgb[1] == NULL) {
			printf("Failed to allocate memory required, aborted.\n");
			if (anim) *gif = NULL;
			return 1;
		}
		rgb[2] = malloc(nb_pixels);
		if (rgb[2] == NULL) {
			printf("Failed to allocate memory required, aborted.\n");
			if (anim) *gif = NULL;
			return 1;
		}
	}

	pente = computePente(&lo, &hi);
	
	for (i = 0; i <= USHRT_MAX; i++) {
		map[i] = round_to_BYTE((float) i * pente);
		if (map[i] == UCHAR_MAX)
			break;
	}
	if (i != USHRT_MAX + 1) {
		/* no more computation needed, just fill with max value */
		for (++i; i <= USHRT_MAX; i++)
			map[i] = UCHAR_MAX;
	}
	
	/* doing the WORD to BYTE conversion, bottom-up */
	if (fit->naxes[2] == 1) {
		int x, y;
		WORD *src = fit->data;
		BYTE *dst = gif_buffer;
		for (y = 0; y < fit->ry; y++) {
			int desty = fit->ry - y - 1;
			int srcpixel = y * fit->rx;
			int dstpixel = desty * fit->rx;
			for (x = 0; x < fit->rx; x++, srcpixel++, dstpixel++) {
				// linear scaling
				tmp_pixel_value = src[srcpixel] - lo;
				dst[dstpixel] = map[tmp_pixel_value];
			}
		}
	} else {
		int channel;
		for (channel = 0; channel < 3; channel++) {
			int x, y;
			WORD *src = fit->pdata[channel];
			BYTE *dst = rgb[channel];
			for (y = 0; y < fit->ry; y++) {
				int desty = fit->ry-y-1;
				int srcpixel = y * fit->rx;
				int dstpixel = desty * fit->rx;
				for (x = 0; x < fit->rx; x++, srcpixel++, dstpixel++) {
					// linear scaling
					if (src[srcpixel] - lo < 0)
						tmp_pixel_value = 0;
					else 	tmp_pixel_value = src[srcpixel] - lo;
					dst[dstpixel] = map[tmp_pixel_value];
				}
			}
		}
	}

	/* transform to GIF data */
	ColorMapObject *OutputColorMap = NULL;
	int ColorMapSize = 256;

	if (fit->naxes[2] == 1) {
		GifColorType gif_colors[256];	// the palette (a.k.a. color map)
		for (i = 0; i < 256; i++) {
			gif_colors[i].Red = gif_colors[i].Green = gif_colors[i].Blue = i;
		}

		if ((OutputColorMap = GifMakeMapObject(ColorMapSize, gif_colors)) == NULL) {
			free(gif_buffer);
			printf("Failed to allocate memory required, aborted.\n");
			if (anim) *gif = NULL;
			return 1;
		}
	} else {

		if ((OutputColorMap = GifMakeMapObject(ColorMapSize, NULL)) == NULL ||
				(gif_buffer = malloc(nb_pixels * sizeof(GifByteType))) == NULL) {
			free(rgb[0]); free(rgb[1]); free(rgb[2]);
			printf("Failed to allocate memory required, aborted.\n");
			if (anim) *gif = NULL;
			return 1;
		}

		if (GifQuantizeBuffer(fit->rx, fit->ry, &ColorMapSize,
					rgb[0], rgb[1], rgb[2],
					gif_buffer, OutputColorMap->Colors) == GIF_ERROR) {
			free(rgb[0]); free(rgb[1]); free(rgb[2]); free(gif_buffer);
			siril_log_message(_("Failed to convert data into GIF data format\n"));
			if (anim) *gif = NULL;
			GifFreeMapObject(OutputColorMap);
			return 1;
		}

		free(rgb[0]); free(rgb[1]); free(rgb[2]);
	}

	/*********** code from SaveGif function **********/
	int error;
	GifFileType *GifFile;
	if (!anim || *gif == NULL) {
		if ((GifFile = EGifOpenFileName(filename, FALSE, &error)) == NULL) {
			siril_log_message(_("Error opening GIF file: %s\n"), GifErrorString(error));
			free(gif_buffer);
			GifFreeMapObject(OutputColorMap);
			return 1;
		}

		if (EGifPutScreenDesc(GifFile, fit->rx, fit->ry, 256, 0, OutputColorMap) == GIF_ERROR) {
			siril_log_message(_("Error describing GIF file: %s\n"), GifErrorString(error));
			EGifCloseFile(GifFile, NULL);
			free(gif_buffer);
			GifFreeMapObject(OutputColorMap);
			return 1;
		}

		if (EGifPutComment(GifFile, "Made with Siril, for free astronomy") == GIF_ERROR) {
			siril_log_message(_("Error adding comment in GIF file\n"));
		}

		/* add loop, optional */
		if (anim && loop_count > 0) {
			char nsle[12] = "NETSCAPE2.0";
			char subblock[3];
			subblock[0] = 1;
			subblock[1] = loop_count / 256;
			subblock[2] = loop_count % 256;
			if (EGifPutExtensionLeader(GifFile, APPLICATION_EXT_FUNC_CODE) == GIF_ERROR ||
					EGifPutExtensionBlock(GifFile, 11, nsle) == GIF_ERROR ||
					EGifPutExtensionBlock(GifFile, 3, subblock) == GIF_ERROR ||
					EGifPutExtensionTrailer(GifFile) == GIF_ERROR) {
				siril_log_message(_("Error adding the loop in GIF file\n"));
			}
		}
	} else {
		GifFile = *gif;
	}

	/* for each image */

	/* add delay, optional */
	if (anim && delay > 0) {
		unsigned char ExtStr[4] = { 0x04, 0x00, 0x00, 0xff };
		ExtStr[1] = delay % 256;
		ExtStr[2] = delay / 256;
		if (EGifPutExtension(GifFile, GRAPHICS_EXT_FUNC_CODE, 4, ExtStr) == GIF_ERROR) {
			siril_log_message(_("Error adding the delay in GIF file\n"));
		}
	}

	if (EGifPutImageDesc(GifFile, 0, 0, fit->rx, fit->ry, FALSE, OutputColorMap) == GIF_ERROR) {
		siril_log_message(_("Error describing GIF file: %s\n"), GifErrorString(error));
		EGifCloseFile(GifFile, NULL);
		free(gif_buffer);
		if (anim) *gif = NULL;
		GifFreeMapObject(OutputColorMap);
		return 1;
	}

	GifByteType *ptr = gif_buffer;
	for (i = 0; i < fit->ry; i++) {
		if (EGifPutLine(GifFile, ptr, fit->rx) == GIF_ERROR) {
			EGifCloseFile(GifFile, NULL);
			free(gif_buffer);
			if (anim) *gif = NULL;
			GifFreeMapObject(OutputColorMap);
			return 1;
		}
		ptr += fit->rx;
	}
	/* end of for each image */

	free(gif_buffer);
	GifFreeMapObject(OutputColorMap);

	if (!anim) {
		if (EGifCloseFile(GifFile, &error) == GIF_ERROR) {
			siril_log_message(_("Error closing GIF file: %s\n"), GifErrorString(error));
			return 1;
		}
	} else {
		*gif = GifFile;
	}
	return 0;
}

void closegif(GifFileType **gif) {
	int error;
	if (gif == NULL || *gif == NULL)
		return;

	if (EGifCloseFile(*gif, &error) == GIF_ERROR) {
		siril_log_message(_("Error closing GIF file: %s\n"), GifErrorString(error));
	}
	*gif = NULL;
}
#endif

