#ifndef SIRIL_H
#define SIRIL_H
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

// PATH_MAX is not available on Hurd at least
#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define SQR(x) ((x)*(x))

#include <gtk/gtk.h>
#include <fitsio.h>	// fitsfile
#include <gsl/gsl_histogram.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define USHRT_MAX_DOUBLE ((double)USHRT_MAX)
#define UCHAR_MAX_DOUBLE ((double)UCHAR_MAX)
#define UCHAR_MAX_SINGLE ((float)UCHAR_MAX)

#define SEQUENCE_DEFAULT_INCLUDE TRUE	// select images by default

typedef unsigned char BYTE;		// default type for image display data
typedef unsigned short WORD;		// default type for internal image data

#define MAX_COMMAND_WORDS 16		// max number of words to split in command line input

#define CMD_HISTORY_SIZE 50		// size of the command line history

#define ZOOM_MAX	16.0
#define ZOOM_MIN	0.0625
#define ZOOM_NONE	1.0
#define ZOOM_FIT	-1.0	// or any value < 0
#define ZOOM_DEFAULT	ZOOM_FIT

/* when requesting an image redraw, it can be asked to remap its data before redrawing it.
 * REMAP_NONE	doesn't remaps the data,
 * REMAP_ONLY	remaps only the current viewport (color channel) and the mixed (RGB) image
 * REMAP_ALL	remaps all view ports, useful when all the colors come from the same file.
 */
#define REMAP_NONE	0
#define REMAP_ONLY	1
#define REMAP_ALL	2

enum {
	COLUMN_FILENAME,		// string
	COLUMN_DATE,		// string
	N_COLUMNS_CONVERT
};


typedef enum {
	TYPEUNDEF,
	TYPEFITS,
	TYPETIFF,
	TYPEBMP,
	TYPEPNG,
	TYPEJPG,
	TYPEPNM,
	TYPEPIC,
	TYPERAW,
	TYPEAVI,
	TYPESER,
} image_type;

#define USE_DARK	0x01
#define USE_FLAT	0x02
#define USE_OFFSET	0x04

/* cookies for the file chooser */
#define OD_NULL 	0
#define OD_FLAT 	1
#define OD_DARK 	2
#define OD_OFFSET 	3
#define OD_CWD 		4
#define OD_OPEN 	5
#define OD_CONVERT 	6

/* indices of the image data layers */
#define BW_LAYER 	0
#define RLAYER		0
#define GLAYER		1
#define BLAYER		2
#define RGB_LAYER 	3

/* indices of the viewports (graphical elements) */
#define BW_VPORT 	0
#define RED_VPORT 	0
#define GREEN_VPORT 	1
#define BLUE_VPORT 	2
#define RGB_VPORT 	3
#define MAXGRAYVPORT 	3			// 3 gray vports supported only (R, G, B)
#define MAXCOLORVPORT	1			// 1 color vport supported only (RGB)
#define MAXVPORT 	MAXGRAYVPORT + MAXCOLORVPORT

/* defines for copyfits actions */
#define CP_INIT		0x01	// initialize data array with 0s
#define CP_ALLOC	0x02	// reallocs data array
#define CP_COPYA	0x04	// copy data array content
#define CP_FORMAT	0x08	// copy size and bitpix info
#define CP_EXTRACT	0x10	// extract a 16bit plane from a 48 bit fit
#define CP_EXPAND	0x20	// expands a 16bit fits to a 48bit one.

#define CONVALL	(1 << 0)	// all files, not based on an index
/* file types available */
#define CONVBMP	(1 << 1)
#define CONVPNG	(1 << 2)
#define CONVJPG	(1 << 3)
#define CONVAVI	(1 << 4)
#define CONVPNM	(1 << 5)
#define CONVRAW (1 << 6)
#define CONVCFA (1 << 7)
#define CONVTIF (1 << 8)
#define CONVPIC (1 << 9)
/* layer conversion type */
#define CONV1X3	(1 << 16)
#define CONV3X1	(1 << 17)
#define CONV1X1	(1 << 18)
#define CONVUFL	(1 << 20)	// use fixed length

/* operations on image data */
#define OPER_ADD 'a'
#define OPER_SUB 's'
#define OPER_MUL 'm'
#define OPER_DIV 'd'

#define PREVIEW_NB 2

#define RESULT_IMAGE -1		// used as current image index in a sequence
				// when the result of a processing is displayed
#define UNRELATED_IMAGE -2	// image loaded while a sequence was loaded too

#define MAX_STARS 20000

/* constants for loglut function */
#define LOG 1
#define EXP -1

#define PROGRESS_NONE -2.0		// don't update the progress bar value
#define PROGRESS_PULSATE -1.0		// pulsate the progress bar
#define PROGRESS_RESET 0.0		// reset the progress bar
#define PROGRESS_DONE 1.0		// fill the progress bar
#define PROGRESS_TEXT_RESET ""		// reset the progress bar's text

typedef struct imdata imgdata;
typedef struct registration_data regdata;
typedef struct layer_info_struct layer_info;
typedef struct sequ sequence;
typedef struct single_image single;
typedef struct ffit fits;
typedef struct libraw_config libraw;
typedef struct stack_config stackconf;
typedef struct cominf cominfo;
typedef struct image_stats imstats;
typedef struct fwhm_struct fitted_PSF;
typedef struct rectangle_struct rectangle;
typedef struct point_struct point;
typedef struct gradient_struct gradient;
typedef struct historic_struct historic;

/* global structures */

/* image data, exists once for each image */
struct imdata {
	int filenum;		/* real file index in the sequence, i.e. for mars9.fit = 9 */
	gboolean incl;		/* selected in the sequence, included for future processings? */
	//double eval;		/* an evaluation of the quality of the image, not yet used */
	/* for deep-sky images, the quality is the FWHM of a star or an average of FWHM of stars,
	 * for planetary and deep-sky images, it's computed from entropy, contrast or other eval functions. */
	imstats *stats;		/* statistics of the image, used as a cache for full image first
				   channel statistics, particularly used in stacking normalization.
				   NULL if not available, this is stored in seqfile if non NULL */
};

/* Procedure signature for sequence processing.
 * Returns < 0 for an error that should stop the processing on the sequence.
 * Other return values are not used.
 * Processed data should be written in the sequence data directly. */
typedef int (*sequence_proc)(sequence *seq, int seq_layer, int frame_no, fits *fit);

/* preprocessing data from GUI */
struct preprocessing_data {
	struct timeval t_start;
	gboolean autolevel;
	gboolean use_ccd_formula;
	float normalisation;
	int retval;
};

/* registration data, exists once for each image and each layer */
struct registration_data {
	int shiftx, shifty;	// we could have a subpixel precision, but is it needed? saved
	float rot_centre_x, rot_centre_y;	// coordinates for the rotation centre, saved
	float angle;		// angle for the rotation, saved
	fitted_PSF *fwhm_data;	// used in PSF/FWHM registration, not saved
	float fwhm;		// copy of fwhm->fwhmx, used as quality indicator, saved data
	//double entropy;		// used in DFT registration, saved data
	// entropy could be replaced by a more general quality indicator at
	// some point. It's the only double value stored in the .seq files, others are single.
	double quality;
};

/* ORDER OF POLYNOMES */
typedef enum {
	POLY_1,
	POLY_2,
	POLY_3,
	POLY_4,
} poly_order;

typedef enum {
	NORMAL_DISPLAY,	
	LOG_DISPLAY,
	SQRT_DISPLAY,
	SQUARED_DISPLAY,
	ASINH_DISPLAY,
	HISTEQ_DISPLAY
} display_mode;			// used in the layer_info_struct below
#define DISPLAY_MODE_MAX HISTEQ_DISPLAY

typedef enum {
	NORMAL_COLOR,	
	RAINBOW_COLOR
} color_map;

typedef enum {
	MIPSLOHI,
	MINMAX,
	USER
} sliders_mode;

typedef enum {
	FILE_CONVERSION,
	IMAGE_SEQ,
	PRE_PROC,
	REGISTRATION,
	STACKING,
	OUTPUT_LOGS
} main_tabs;

typedef enum {
	BAYER_BILINEAR,
	BAYER_NEARESNEIGHBOR,
	BAYER_VNG,
	BAYER_AHD,
	BAYER_SUPER_PIXEL,
} interpolation_method;

typedef enum {
    BAYER_FILTER_RGGB,
    BAYER_FILTER_BGGR,
    BAYER_FILTER_GBRG,
    BAYER_FILTER_GRBG,
    BAYER_FILTER_NONE = -1		//case where bayer pattern is undefined or untested
} sensor_pattern ;
#define BAYER_FILTER_MIN BAYER_FILTER_RGGB
#define BAYER_FILTER_MAX BAYER_FILTER_GRBG

struct layer_info_struct {
	char *name;			// name of the layer (a filter name)
	double wavelength;		// the wavelength of the filter, in nanometres
	WORD lo, hi;			// the values of the cutoff sliders
	//WORD min, max;		// the min and max values of the sliders
	gboolean cut_over, cut_under;	// display values over hi or under lo as negative
	display_mode rendering_mode;	// defaults to NORMAL_DISPLAY
};

typedef enum { SEQ_REGULAR, SEQ_SER,
#ifdef HAVE_FFMS2
	SEQ_AVI,
#endif
	SEQ_INTERNAL
} sequence_type;

struct sequ {
	char *seqname;		// name of the sequence, as in name.seq
	int number;		// number of images in the sequence
	int selnum;		// number of selected images in the sequence
	int fixed;		// fixed length of image index in filename (like %3d)
	int nb_layers;		// number of layers embedded in each image file
	unsigned int rx;	// first image width
	unsigned int ry;	// first image height
	layer_info *layers;	// info about layers
	int reference_image;	// reference image for registration
	imgdata *imgparam;	// a structure for each image of the sequence
	regdata **regparam;	// *regparam[nb_layers]
	/* beg and end are used prior to imgparam allocation, hence their usefulness */
	int beg;		// imgparam[0]->filenum
	int end;		// imgparam[number-1]->filenum

	/* Data used when a new sequence is builded
	 * --> in global star registration for example
	 */
	int new_total;

	/* registration previsualisation and manual alignment data */
	int previewX[PREVIEW_NB], previewY[PREVIEW_NB];	// center, -1 is uninitialized value
	int previewW[PREVIEW_NB], previewH[PREVIEW_NB];	// 0 is uninitialized value

	sequence_type type;
	struct ser_struct *ser_file;
#ifdef HAVE_FFMS2
	struct film_struct *film_file;
	char *ext;		// extension of video, NULL if not video
#endif
	fits **internal_fits;	// for INTERNAL sequences: images references. Length: number
	fitsfile **fptr;	// file descriptors for open-mode operations
#ifdef _OPENMP
	omp_lock_t *fd_lock;	// locks for open-mode threaded operations
#endif

	fits *offset;		// the image containing offset data
	fits *dark;		// the image containing dark data
	fits *flat;		// the image containing flat data
	char *ppprefix;		// prefix for filename output of preprocessing
	int current;		// file number currently loaded in wfit (or displayed)
	//struct registration_method reg_method;	// is it the right place for that?
	
	gboolean needs_saving;	// a dirty flag for the sequence, avoid saving it too often
};

/* this struct is used to manage data associated with a single image loaded, outside a sequence */
struct single_image {
	char *filename;		// the name of the file
	char *comment;		// comment on how the file got there (user load, result...)
	int nb_layers;		// number of layers embedded in each image file
	layer_info *layers;	// info about layers
	fits *fit;		// the fits is still gfit, but a reference doesn't hurt

	/* enabling pre-processing on a single image */
	fits *offset;		// the image containing offset data
	fits *dark;		// the image containing dark data
	fits *flat;		// the image containing flat data
	char *ppprefix;		// prefix for filename output of preprocessing
};

struct ffit {
	unsigned int rx;	// image width	(naxes[0])
	unsigned int ry;	// image height	(naxes[1])
	int bitpix;
	/* bitpix can take the following values:
	 * BYTE_IMG	(8-bit byte pixels, 0 - 255)
	 * SHORT_IMG	(16 bit signed integer pixels)	
	 * USHORT_IMG	(16 bit unsigned integer pixels)	(used by Siril, quite off-standard)
	 * LONG_IMG	(32-bit integer pixels)
	 * FLOAT_IMG	(32-bit floating point pixels)
	 * DOUBLE_IMG	(64-bit floating point pixels)
	 * http://heasarc.nasa.gov/docs/software/fitsio/quick/node9.html
	 */
	int naxis;		// number of dimensions of the image
	long naxes[3];		// size of each dimension
	/* naxes[0] is rx, naxes[1] is ry
	 * Then, for gray images, naxes[2] is unused in FITS but set to 1, and naxes is 2.
	 * For RGB images, naxes[2] is 3 and naxis is 3.
	 * */

	/* data obtained from the FITS file */
	WORD lo;	// MIPS-LO key in FITS file, which is "Lower visualization cutoff"
	WORD hi;	// MIPS-HI key in FITS file, which is "Upper visualization cutoff"
	float pixel_size_x, pixel_size_y;	// XPIXSZ and YPIXSZ keys
	short binning_x, binning_y;		// XBINNING and YBINNING keys
	char date_obs[FLEN_VALUE];		// YYYY-MM-DDThh:mm:ss observation start, UT
	char instrume[FLEN_VALUE];		// INSTRUME key
	/* data obtained from FITS or RAW files */
	double focal_length, iso_speed, exposure, aperture, ccd_temp;

	/* data used in the Fourrier space */
	double dft_norm[3];			// Normalization value
	char dft_type[FLEN_VALUE];		// spectrum, phase
	char dft_ord[FLEN_VALUE];		// regular, centered
	unsigned int dft_rx, dft_ry;		// padding: original value of picture size
	
	/* data computed or set by Siril */
	unsigned short min[3];	// min for all layers
	unsigned short max[3];	// max for all layers
	unsigned short maxi;	// max of the max[3]
	unsigned short mini;	// min of the min[3]

	fitsfile *fptr;		// file descriptor. Only used for file read and write.
	WORD *data;		// 16-bit image data (depending on image type)
	WORD *pdata[3];		// pointers on data, per layer data access (RGB)
	char *header;		// entire header of the FITS file. NULL for non-FITS file.
};

/* This structure is used for all the elements in the box libraw_settings.
 * Don't forget to update conversion.c:initialize_libraw_settings() data when
 * modifying the glade settings */
struct libraw_config {
	gboolean cfa, ser_cfa, ser_force_bayer;
	double mul[3], bright;					// Color  & brightness adjustement mul[0] = red, mul[1] = green = 1, mul[2] = blue
	int auto_mul, use_camera_wb, use_auto_wb;		// White Balance parameters
	int user_qual;						// Index of the Matrix interpolation set in dcraw, 0: bilinear, 1: VNG, 2: PPG, 3: AHD
	int user_black;						// black point correction
	double gamm[2];						// Gamma correction
	sensor_pattern bayer_pattern;
	interpolation_method bayer_inter;					
};

struct stack_config {
	int method;				// 0=sum, 1=median, 2=average, 3=pixel max - Use to save preferences in the init file
	int normalisation_method;
	int rej_method;
	double memory_percent;			// percent of available memory to use for stacking
};

struct rectangle_struct {
	int x, y, w, h;
};

struct point_struct {
	double x, y;
};

struct gradient_struct {
	point centre;
	double boxvalue;
};

struct historic_struct {
	char *filename;
	char history[FLEN_VALUE];
	int rx, ry;
};

struct cominf {
	/* current version of GTK, through GdkPixmap, doesn't handle gray images, so
	 * graybufs are the same size than the rgbbuf with 3 times the same value */
	guchar *graybuf[MAXGRAYVPORT];	// one B/W display buffer per viewport (R,G,B)
	guchar *rgbbuf;			// one rgb display buffer
	/* cairo image surface related data */
	int surface_stride[MAXVPORT];	// allocated stride
	int surface_height[MAXVPORT];	// allocated height
	cairo_surface_t *surface[MAXVPORT];
	gboolean buf_is_dirty[MAXVPORT];// dirtyness of each buffer (= need to redraw)
	
	/* Color map */
	color_map color;

	GtkWidget *vport[MAXVPORT];	// one drawingarea per layer, one rgb drawingarea
	int cvport;			// current viewport, index in the list vport above
	GtkAdjustment *hadj[MAXVPORT];	// adjustments of vport scrollbars
	GtkAdjustment *vadj[MAXVPORT];	// adjustments of vport scrollbars
	sliders_mode sliders;		// 0: min/max, 1: MIPS-LO/HI, 2: user
	gboolean leveldrag;		// middle click being dragged if true
	int preprostatus;
	int preproformula;		// 0: APN, 1: CCD
	gboolean show_excluded;		// show excluded images in sequences NOT USED!
	double zoom_value;		// 1.0 is normal zoom, use get_zoom_val() to access it

	/* selection rectangle for registration, FWHM, PSF */
	gboolean drawn;			// true if the selection rectangle has been drawn TO REMOVE!
	gboolean drawing;		// true if the rectangle is being set (clicked motion)
	gint startX, startY;		// where the mouse was originally clicked to
	rectangle selection;		// coordinates of the selection rectangle

	/* alignment preview data */
	//guchar *preview_buf[PREVIEW_NB];
	cairo_surface_t *preview_surface[PREVIEW_NB];
	GtkWidget *preview_area[PREVIEW_NB];
	guchar *refimage_regbuffer;	// the graybuf[registration_layer] of the reference image
	cairo_surface_t *refimage_surface;

	char *wd;			// working directory, where images and sequences are
	char *initfile;			// the path of the init file
	
	char *ext;		// FITS extension used in SIRIL
	int len_ext;

	int reg_settings;		// 0=DFT Translation, 1=PSF Translation - Use to save preferences in the init file
	
	stackconf stack;
	
	int filter;

	/* history of the command line. This is a circular buffer (cmd_history)
	 * of size cmd_hist_size, position to be written is cmd_hist_current and
	 * position being browser for display of the history is cmd_hist_display.
	 */
	char **cmd_history;		// the history of the command line
	int cmd_hist_size;		// allocated size
	int cmd_hist_current;		// current command index
	int cmd_hist_display;		// displayed command index

	/* history of operations */
	historic *history;			// the history of all operations
	int hist_size;			// allocated size
	int hist_current;		// current index
	int hist_display;		// displayed index
	char *swap_dir;

	libraw raw_set;			// the libraw setting

	sequence seq;			// currently loaded sequence	TODO: *seq
	single *uniq;			// currently loaded image, if outside sequence

	gsl_histogram *layers_hist[MAXVPORT]; // current image's histograms

	fitted_PSF **stars;		// list of stars detected in the current image
	gboolean star_is_seqdata;	// the only star in stars belongs to seq, don't free it
	int selected_star;		// current selected star in the GtkListStore
	
	gradient *grad;
	int grad_nb_boxes, grad_size_boxes;
	gboolean grad_boxes_drawn;

	GThread *thread;		// the thread for processings
	GMutex mutex;			// a mutex we use for this thread
	gboolean run_thread;		// the main thread loop condition
	int max_thread;		// maximum of thread used
};

/* this structure is used to characterize the statistics of the image */
struct image_stats {
	double mean, avgdev, median, sigma, min, max;
	char layername[6];
};

#if 0
/* TODO: this structure aims to allow the composition of several 1-channel images and make
 * more easy the management of RGB compositing */
typedef struct image_layer_struct image_layer;
struct image_layer_struct {
	char		*layer_name;		/* the name of the layer (the color or band name) */
	fits		*fit;			/* fits data of the layer */
	int		naxis;			/* this image is naxis in fits */
	guchar		*graybuf;		/* mapped image for display purposes */
	int		stride;			/* Cairo data width */
	cairo_surface_t	*surface;		/* Cairo image surface */
	GtkWidget	*vport;			/* the viewport */
	double		wavelength;		/* the wavelength associated with the channel */
	guchar		rmap, gmap, bmap;	/* mapping to rgb colors, initialized function of the wavelength */
	unsigned int	hi, lo;			/* same as fits_data->{hi,lo} but for display/compositing purposes */
	char		*filename;		/* the filename of the fits_data file */
	int		naxis;			/* axis number of the fits file filename, may not be 0 if it's an RGB fits for example */
};
#endif

#ifndef MAIN
extern GtkBuilder *builder;	// get widget references anywhere
extern cominfo com;		// the main data struct
extern fits gfit;		// currently loaded image
extern fits wfit[5];		// used for temp files, can probably be replaced by local variables
extern char **supported_extensions;
#endif

#endif /*SIRIL */