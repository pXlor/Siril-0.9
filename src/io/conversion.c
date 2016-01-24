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
*/

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <gdk/gdkkeysyms.h>
#include <gtk/gtk.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <dirent.h>
#include "core/siril.h"
#include "core/proto.h"
#include "io/conversion.h"
#include "io/films.h"
#include "gui/callbacks.h"
#include "algos/demosaicing.h"

#define SUFFIX_MAX_LEN 9
#define MAX_OF_EXTENSIONS 50

//~ static char *sourceroot = NULL;
static char *destroot = NULL;
static char sourcesuf[SUFFIX_MAX_LEN+1];	// initialized at runtime
static char destsuf[]="fit";
static unsigned int convflags = CONV1X3 | CONVBMP | CONVPIC | CONVTIF | CONVJPG | CONVPNG | CONVRAW | CONVCFA | CONVALL;
static unsigned int supported_filetypes = 0;	// initialized by initialize_converters()
//static char pnmpath[255]={"/usr/bin"};
//static char pnmcom[3][255]={"pngtopnm","jpegtopnm","bmptopnm"};
static char combuffer[512];

//~ static gpointer convert_thread_worker(gpointer p);
//~ static gboolean end_convert_idle(gpointer p);
static gpointer convert_thread_worker(gpointer p);
static gboolean end_convert_idle(gpointer p);

char **supported_extensions;		// initialized by initialize_converters() and is a NULL-terminated array
supported_raw_list supported_raw[] = {
	
	{"dng",	"Adobe", BAYER_FILTER_RGGB},
	{"mos",	"Aptus", BAYER_FILTER_RGGB},
	{"cr2",	"Canon", BAYER_FILTER_RGGB},
	{"crw",	"Canon", BAYER_FILTER_RGGB},
	{"bay",	"Casio", BAYER_FILTER_NONE},		// Not tested
	{"erf",	"Epson", BAYER_FILTER_RGGB},
	{"raf",	"Fuji", BAYER_FILTER_RGGB},		// Bugged with some files
	{"3fr",	"Hasselblad", BAYER_FILTER_GRBG},	// GRBG, RGGB		
	{"kdc",	"Kodak", BAYER_FILTER_GRBG},
	{"dcr",	"Kodak", BAYER_FILTER_GRBG},
	{"mef",	"Mamiya", BAYER_FILTER_RGGB},
	{"mrw",	"Minolta", BAYER_FILTER_RGGB},
	{"nef",	"Nikon", BAYER_FILTER_RGGB},
	{"nrw",	"Nikon", BAYER_FILTER_RGGB},
	{"orf",	"Olympus", BAYER_FILTER_GRBG},
	{"raw",	"Leica", BAYER_FILTER_RGGB},
	{"rw2",	"Panasonic", BAYER_FILTER_BGGR},
	{"pef",	"Pentax", BAYER_FILTER_BGGR},
	{"ptx",	"Pentax", BAYER_FILTER_NONE},		// Not tested
	{"x3f",	"Sigma", BAYER_FILTER_NONE},		// No Bayer-Pattern
	{"srw",	"Samsung", BAYER_FILTER_BGGR},
	{"arw",	"Sony", BAYER_FILTER_RGGB}
};

int get_nb_raw_supported() {
	return sizeof(supported_raw) / sizeof(supported_raw_list);
}
		
char *conversion_tips[] = {
	//~ Color interpolation
	"Your RAW files are being converted to interpolated colour-pictures. "
	"It would be wrong to use this mode for pre-processing. "
	"To change this behaviour, check the correspondant button in File->Settings->Raw images.",
	//~ CFA conversion
	"Your RAW files are being converted to CFA monochrome pictures. "
	"To change this behaviour, uncheck the correspondant button in File->Settings->Raw images." 
};

char *filter_pattern[] = {
	"RGGB",
	"BGGR",
	"GBRG",
	"GRBG"
};

	/* This function is used with command line only */ 
void list_format_available() {
	puts("======================================================="); 
	puts("[            Supported image file formats             ]");
	puts("======================================================="); 
	puts("FITS\t(*.fit, *.fits, *.fts)");
	puts("BMP\t(*.bmp)");
	puts("NetPBM\t(*.ppm, *.pgm, *.pnm)");
	puts("PIC\t(*.pic)");
#ifdef HAVE_LIBRAW
	printf("RAW\t(");
	int i, nb_raw;
	
	nb_raw = get_nb_raw_supported();
	for (i = 0; i < nb_raw; i++) {
		printf("*.%s",supported_raw[i].extension);
		if (i != nb_raw - 1) printf(", ");
	}
	printf(")\n");
#endif

#ifdef HAVE_LIBTIFF
	puts("TIFF\t(*.tif, *.tiff)");
#endif
#ifdef HAVE_LIBJPEG
	puts("JPEG\t(*.jpg, *.jpeg)");
#endif
#ifdef HAVE_LIBPNG
	puts("PNG\t(*.png)");
#endif
}

void check_for_conversion_form_completeness() {
	static GtkTreeView *tree_convert = NULL;
	GtkTreeIter iter;
	GtkTreeModel *model = NULL;
	gboolean valid;
	GtkWidget *go_button = lookup_widget("convert_button");
	
	if (tree_convert == NULL)
		tree_convert = GTK_TREE_VIEW(gtk_builder_get_object(builder, "treeview_convert"));
	
	model = gtk_tree_view_get_model(tree_convert);
	valid = gtk_tree_model_get_iter_first(model, &iter);
	gtk_widget_set_sensitive (go_button, destroot && destroot[0] != '\0' && valid);
	update_statusbar_convert();
}

void on_convtoroot_changed (GtkEditable *editable, gpointer user_data){
	const char *name=gtk_entry_get_text(GTK_ENTRY(editable));
	if (destroot) free(destroot);
	destroot = strdup(name);
	check_for_conversion_form_completeness();
}

void on_conv3planefit_toggled (GtkToggleButton *togglebutton, gpointer user_data){
	convflags |= CONV1X3;
	convflags &= ~(CONV3X1|CONV1X1);
}

void on_conv3_1plane_toggled (GtkToggleButton *togglebutton, gpointer user_data) {
	convflags |= CONV3X1;
	convflags &= ~(CONV1X1|CONV1X3);
}

void on_radiobutton_conv_cfa_toggled (GtkToggleButton *togglebutton, gpointer user_data) {
	static GtkTreeView *tree_convert = NULL;
	GtkTreeIter iter;
	GtkTreeModel *model = NULL;
	gboolean valid;
	
	if (tree_convert == NULL)
		tree_convert = GTK_TREE_VIEW(gtk_builder_get_object(builder, "treeview_convert"));
	model = gtk_tree_view_get_model(tree_convert);
	valid = gtk_tree_model_get_iter_first(model, &iter);
	com.is_cfa = gtk_toggle_button_get_active(togglebutton);
	
	while(valid) {
		gchar *str_data, *msg;
		gtk_tree_model_get (model, &iter, 0, &str_data, -1);	//0 for FILECOLUMN
		GtkToggleButton *but = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_conv"));
		if (!ends_with(str_data, ".fit") && (com.is_cfa)) {
			g_signal_handlers_block_by_func(but, on_radiobutton_conv_cfa_toggled, NULL);
			gtk_toggle_button_set_active (but, TRUE);
			g_signal_handlers_unblock_by_func(but, on_radiobutton_conv_cfa_toggled, NULL);
			msg = siril_log_message("There is at least one file not compatible. For demosaicing, all files must have *.fit extension\n");
			show_dialog(msg, "ERROR", "gtk-dialog-error");
			return;
		}
		/*else if (ends_with(str_data, ".fit") && (!is_cfa)) {
			g_signal_handlers_block_by_func(togglebutton, on_radiobutton_conv_cfa_toggled, NULL);
			gtk_toggle_button_set_active (togglebutton, TRUE);
			g_signal_handlers_unblock_by_func(togglebutton, on_radiobutton_conv_cfa_toggled, NULL);
			msg = siril_log_message("There is at least one file not compatible. No *.fit files allowed in FITS conversion.\n");
			show_dialog(msg, "ERROR", "gtk-dialog-error");
			return;
		}*/
		valid = gtk_tree_model_iter_next (model, &iter);
	}
}

/* COMMAND LINE CONVERSION UTILITIES */

/* returns true if the command mplayer for AVI file importing is available */
gboolean mplayer_is_available() {
	int retval = system("mplayer > /dev/null 2>&1");
	if (WIFEXITED(retval))
		return 0 == WEXITSTATUS(retval); // mplayer always returns 0
	return FALSE;
}

/*************************************************************************************
 * 
 * 
 * ********************************************************************************/
 
void update_raw_cfa_tooltip() {
	GtkWidget* conv_button = lookup_widget("convert_button");
	gtk_widget_set_tooltip_text(conv_button, conversion_tips[(int)com.raw_set.cfa]);
}

/* This function sets all default values of libraw settings in the com.raw_set
 * struct, as defined in the glade file.
 * When the ini file is read, the values of com.raw_set are overwritten, but if the
 * file is missing, like the first time Siril is launched, we don't want to have the
 * GUI states reset to zero by set_GUI_LIBRAW() because the data in com.raw_set had
 * not been initialized with the default GUI values (= initialized to 0).
 */
void initialize_libraw_settings() {
	com.raw_set.cfa = TRUE;		// CFA
	com.raw_set.bright = 1.0;		// brightness
	com.raw_set.mul[0] = 1.0;		// multipliers: red
	com.raw_set.mul[1] = 1.0;		// multipliers: green, not used because always equal to 1
	com.raw_set.mul[2] = 1.0;		// multipliers: blue
	com.raw_set.auto_mul = 1;		// multipliers are Either read from file, or calculated on the basis of file data, or taken from hardcoded constants
	com.raw_set.user_black = 0;		// black point correction
	com.raw_set.use_camera_wb = 0;	// if possible, use the white balance from the camera. 
	com.raw_set.use_auto_wb = 0;		// use automatic white balance obtained after averaging over the entire image
	com.raw_set.user_qual = 1;		// type of interpolation. AHD by default
	com.raw_set.gamm[0] = 1.0;		// gamm curve: linear by default
	com.raw_set.gamm[1] = 1.0;
	com.raw_set.bayer_pattern = BAYER_FILTER_RGGB;
	com.raw_set.bayer_inter = BAYER_VNG;
}

void initialize_ser_debayer_settings() {
	com.raw_set.ser_cfa = TRUE;
	com.raw_set.ser_force_bayer = TRUE;
}
 
/* initialize converters (utilities used for different image types importing) *
 * updates the label listing the supported input file formats, and modifies the
 * list of file types used in convflags */
void initialize_converters() {
	GtkLabel *label_supported;
	gchar text[256];
	gboolean has_entry = FALSE;	// true if something is supported
	sprintf(text, "\t");		// if changed, change its removal before log_message
	//text[0] = '\0';
	int count_ext = 0;
	/* internal converters */
	supported_filetypes |= CONVBMP;
	strcat(text, "BMP images, ");
	supported_filetypes |= CONVPIC;
	strcat(text, "PIC images (IRIS), ");
	supported_filetypes |= CONVPNM;
	strcat(text, "PGM and PPM binary images");
	has_entry = TRUE;
		
	supported_extensions = malloc(MAX_OF_EXTENSIONS * sizeof(char*));
	/* internal extensions */
	if (supported_extensions == NULL) {
		fprintf(stderr, "initialize_converters: error allocating data\n");
		return;
	}
	supported_extensions[count_ext++] = ".fit";
	supported_extensions[count_ext++] = ".fits";
	supported_extensions[count_ext++] = ".fts";
	supported_extensions[count_ext++] = ".bmp";
	//~ supported_extensions[count_ext++] = ".ser";			// it is like avi, it is not use here for now
	supported_extensions[count_ext++] = ".ppm";
	supported_extensions[count_ext++] = ".pgm";
	supported_extensions[count_ext++] = ".pnm";
	supported_extensions[count_ext++] = ".pic";
	
	initialize_ser_debayer_settings();	// below in the file

#ifdef HAVE_LIBRAW
	int i, nb_raw;
	
	supported_filetypes |= CONVRAW;
	if (has_entry)	strcat(text, ", ");
	strcat(text, "RAW images");
	has_entry = TRUE;
	set_libraw_settings_menu_available(TRUE);	// enable libraw settings
	initialize_libraw_settings();	// below in the file
	
	nb_raw = get_nb_raw_supported();
	for (i = 0; i < nb_raw; i++) {
		supported_extensions[count_ext+i] = malloc(strlen(supported_raw[i].extension) + 2 * sizeof (char));
		strcpy(supported_extensions[count_ext+i], ".");
		strcat(supported_extensions[count_ext+i], supported_raw[i].extension);
	}
	count_ext += nb_raw;
#else
	set_libraw_settings_menu_available(FALSE);	// disable libraw settings
#endif
	supported_filetypes |= CONVCFA;
	if (has_entry)	strcat(text, ", ");
	strcat(text, "FITS-CFA images");
	has_entry = TRUE;

	if (mplayer_is_available() && (supported_filetypes & CONVPNM)) {
		// pnmtofits is used after mplayer extracted frames to PNM
		supported_filetypes |= CONVAVI;
		if (has_entry)	strcat(text, ", ");
		strcat(text, "Videos");
		has_entry = TRUE;
	}

	/* library converters (detected by configure) */
#ifdef HAVE_LIBTIFF
	supported_filetypes |= CONVTIF;
	if (has_entry)	strcat(text, ", ");
	strcat(text, "TIFF images");
	supported_extensions[count_ext++] = ".tif";
	supported_extensions[count_ext++] = ".tiff";
	has_entry = TRUE;
#endif
#ifdef HAVE_LIBJPEG
	supported_filetypes |= CONVJPG;
	if (has_entry)	strcat(text, ", ");
	strcat(text, "JPG images");
	supported_extensions[count_ext++] = ".jpg";
	supported_extensions[count_ext++] = ".jpeg";
	has_entry = TRUE;
#endif
#ifdef HAVE_LIBPNG
	supported_filetypes |= CONVPNG;
	if (has_entry)	strcat(text, ", ");
	strcat(text, "PNG images");
	supported_extensions[count_ext++] = ".png";
	has_entry = TRUE;
#endif
	supported_extensions[count_ext++] = NULL;		// NULL-terminated array

	if (!has_entry)
		sprintf(text, "ERROR: no input file type supported for conversion.");
	else strcat(text, ".");
	label_supported = GTK_LABEL(gtk_builder_get_object(builder, "label9"));
	gtk_label_set_text(label_supported, text);
	siril_log_message("Supported file types: %s\n", text+1);
}

/**************************************************************************/

int check_for_raw_extensions() {
	int i, nb_raw;
	
	nb_raw = get_nb_raw_supported();
	for (i = 0; i < nb_raw; i++) {
		if (!strcasecmp(sourcesuf, supported_raw[i].extension)) return 0;
	}
	return 1;
}

image_type get_type_for_extension_name(const char *extension) {
	if (!strcasecmp(extension, ".fit"))	return TYPEFITS;
	if (!strcasecmp(extension, ".fits"))	return TYPEFITS;
	if (!strcasecmp(extension, ".fts"))	return TYPEFITS;
	if (!strcasecmp(extension, ".png"))	return TYPEPNG;
	if (!strcasecmp(extension, ".jpg"))	return TYPEJPG;
	if (!strcasecmp(extension, ".jpeg"))	return TYPEJPG;
	if (!strcasecmp(extension, ".bmp"))	return TYPEBMP;
	if (!strcasecmp(extension, ".ser"))	return TYPESER;
	if (!strcasecmp(extension, ".tif"))	return TYPETIFF;
	if (!strcasecmp(extension, ".tiff"))	return TYPETIFF;
	if (!strcasecmp(extension, ".pnm"))	return TYPEPNM;
	if (!strcasecmp(extension, ".pgm"))	return TYPEPNM;
	if (!strcasecmp(extension, ".ppm"))	return TYPEPNM;
	if (!strcasecmp(extension, ".pic"))	return TYPEPIC;

	/* Check for raw extentions
	 * Cannot use check_for_raw_extensions() because of the dot */
	int i, nb_raw;
	nb_raw = get_nb_raw_supported();
	for (i = 0; i < nb_raw; i++) {
		char temp[10];
		strncpy(temp, ".", 10);
		strcat(temp, supported_raw[i].extension);
		if (!strcasecmp(extension, temp)) return TYPERAW;
	}

#ifdef HAVE_FFMS2	
	/* Check for film extensions */
	int nb_film;
	nb_film = get_nb_film_ext_supported();
	for (i = 0; i < nb_film; i++) {
		char temp[10];
		strncpy(temp, ".", 10);
		strcat(temp, supported_film[i].extension);
		if (!strcasecmp(extension, temp)) return TYPEAVI;
	}
#else
	if (!strcasecmp(extension, ".avi"))	return TYPEAVI;
#endif
	return TYPEUNDEF;
}

/* takes the extension given by the user and converts it to a convflags flag *
 * returns: 0 is filetype is supported, 1 else */
int set_convflags_from_extension() {
	convflags &= (CONVALL | CONV1X3 | CONV3X1 | CONV1X1 | CONVUFL);	// reset
	if (!strcasecmp(sourcesuf, "bmp")) {
		if (supported_filetypes & CONVBMP) {
			convflags |= CONVBMP;
			return 0;
		}
	} else if (!strcasecmp(sourcesuf, "jpg") || !strcasecmp(sourcesuf, "jpeg")) {
		if (supported_filetypes & CONVJPG) {
			convflags |= CONVJPG;
			return 0;
		}
	} else if (!strcasecmp(sourcesuf, "tif") || !strcasecmp(sourcesuf, "tiff")) {
		if (supported_filetypes & CONVTIF) {
			convflags |= CONVTIF;
			return 0;
		}
	} else if (!strcasecmp(sourcesuf, "png")) { 
		if (supported_filetypes & CONVPNG) {
			convflags |= CONVPNG;
			return 0;
		}
	} else if (!strcasecmp(sourcesuf, "avi") || !strcasecmp(sourcesuf, "mpg") ||
				!strcasecmp(sourcesuf, "mpeg")) {
		if (supported_filetypes & CONVAVI) {	// dependencies already managed
			convflags |= CONVAVI;
			return 0;
		}
	} else if (!strcasecmp(sourcesuf, "pnm") || !strcasecmp(sourcesuf, "ppm") ||
			!strcasecmp(sourcesuf, "pgm")) {
		if (supported_filetypes & CONVPNM) {
			convflags |= CONVPNM;
			return 0;
		}
	} else if (!strcasecmp(sourcesuf, "pic")){
		if (supported_filetypes & CONVPIC) {
			convflags |= CONVPIC;
			return 0;
		}
	} else if (!check_for_raw_extensions()) {
		if (supported_filetypes & CONVRAW) {
			convflags |= CONVRAW;
			return 0;
		}
	} else if (!strcasecmp(sourcesuf, "fit") || !strcasecmp(sourcesuf, "fits") ||
			!strcasecmp(sourcesuf, "fts")) {
		if (com.is_cfa && (supported_filetypes & CONVCFA))
			convflags |= CONVCFA;
		else convflags |= CONVFIT;
		return 0;
	}
	return 1;	// not recognized of not in supported list
}

// appends a '_' to destroot if it ends with a digit
void verify_destroot_validity() {
	if (!destroot) return;
	int len = strlen(destroot);
	if (destroot[len-1] >= '0' && destroot[len-1] <= '9') {
		char *newdest = realloc(destroot, strlen(destroot)+2);
		if (newdest == NULL) {
			free(destroot);
			destroot = NULL;
		} else {
			destroot = newdest;
			destroot[len] = '_';
			destroot[len+1] = '\0';
		}
	}
}

int count_selected_files() {
	static GtkTreeView *tree_convert = NULL;
	GtkTreeModel *model = NULL;
	GtkTreeIter iter;
	gboolean valid;
	int count = 0;
	
	if (tree_convert == NULL)
		tree_convert = GTK_TREE_VIEW(gtk_builder_get_object(builder, "treeview_convert"));

	model = gtk_tree_view_get_model(tree_convert);
	valid = gtk_tree_model_get_iter_first(model, &iter);
	
	while (valid) {
		gchar *file_name, *file_date;
		gtk_tree_model_get (model, &iter, COLUMN_FILENAME, &file_name,
											COLUMN_DATE, &file_date,
											-1);
		valid = gtk_tree_model_iter_next (model, &iter);
		count ++;
	}
	return count;
}

struct _convert_data {
	struct timeval t_start;
	DIR *dir;
	GList *list;
	int start;
	int total;
};

void on_convert_button_clicked(GtkButton *button, gpointer user_data) {
	DIR *dir;
	gchar *file_data, *file_date;
	const gchar *indice;
	static GtkTreeView *tree_convert = NULL;
	static GtkEntry *startEntry = NULL;
	GtkTreeModel *model = NULL;
	GtkTreeIter iter;
	gboolean valid;
	GList *list = NULL;
	int count = 0;
	
	if (tree_convert == NULL) {
		tree_convert = GTK_TREE_VIEW(gtk_builder_get_object(builder, "treeview_convert"));
		startEntry = GTK_ENTRY(gtk_builder_get_object(builder, "startIndiceEntry"));
	}

	struct timeval t_start, t_end;
	
	if (((convflags & CONVALL) || (convflags & CONVAVI)) && get_thread_run()) {
		siril_log_message("Another task is already in progress, ignoring new request.\n");
		return;
	}

	model = gtk_tree_view_get_model(tree_convert);
	valid = gtk_tree_model_get_iter_first(model, &iter);
	if (valid == FALSE) return;	//The tree is empty
	
	while (valid) {
		gtk_tree_model_get (model, &iter, COLUMN_FILENAME, &file_data,
				COLUMN_DATE, &file_date,
				-1);
		list = g_list_append (list, file_data);
		valid = gtk_tree_model_iter_next (model, &iter);
		count ++;
	}
	
	indice = gtk_entry_get_text(startEntry);

	siril_log_color_message("Conversion: processing...\n", "red");	
	gettimeofday(&t_start, NULL);
	
	set_cursor_waiting(TRUE);
	control_window_switch_to_tab(OUTPUT_LOGS);
	
	verify_destroot_validity();
	
	/* then, convert files to Siril's FITS format */

	if (convflags & CONVALL) {
		struct _convert_data *args;
		set_cursor_waiting(TRUE);
		char *tmpmsg;
		if (!com.wd) {
			tmpmsg = siril_log_message("Conversion: no working directory set.\n");
			show_dialog(tmpmsg, "Warning", "gtk-dialog-warning");
			set_cursor_waiting(FALSE);
			return;
		}
		if((dir = opendir(com.wd)) == NULL){
			tmpmsg = siril_log_message("Conversion: error opening working directory %s.\n", com.wd);
			show_dialog(tmpmsg, "Error", "gtk-dialog-error");
			set_cursor_waiting(FALSE);
			return ;
		}

		args = malloc(sizeof(struct _convert_data));
		args->start = (atof(indice) == 0 || atof(indice) > USHRT_MAX) ? 1 : atof(indice);
		args->dir = dir;
		args->list = list;
		args->total = count;
		args->t_start.tv_sec = t_start.tv_sec;
		args->t_start.tv_usec = t_start.tv_usec;
		start_in_new_thread(convert_thread_worker, args);
		return;
	}

	// non-threaded end
	update_used_memory();
	set_cursor_waiting(FALSE);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
}

int convert_film(char *filename) {
	int retval;
	
	set_progress_bar_data("Converting AVI file, it may take some time...", PROGRESS_NONE);
	siril_log_message("Converting AVI file, it may take some time...\n");
	set_progress_bar_data(NULL, 0.05);
	if (convflags & CONV1X1) {
		snprintf(combuffer, 511, "/usr/bin/mplayer -fps 500 -vo pnm:pgm:outdir=%s %s",
				destroot, replace_spaces_from_filename(filename));
		sprintf(sourcesuf, "pgm");
	} else {
		snprintf(combuffer, 511, "/usr/bin/mplayer -fps 500 -vo pnm:outdir=%s %s",
				destroot, replace_spaces_from_filename(filename));
		sprintf(sourcesuf, "ppm");
	}
	fprintf(stdout, "converting avi file using this command: %s\n", combuffer);
	retval = system(combuffer);
	if (retval) return retval;
	if (!WIFEXITED(retval) && 0 != WEXITSTATUS(retval)) {	// mplayer always returns 0
		char *msg = siril_log_message("Invoking mplayer failed, verify you installed it\n");
		set_progress_bar_data(msg, PROGRESS_DONE);
		return -1;
	}
	siril_log_message("Conversion finished, now converting temporary PGM files\n");
	set_progress_bar_data(NULL, 0.35);
	convflags |= (CONVALL | CONVPNM);
	return retval;
}

// This idle function was not really required, but allows to properly join the thread.
static gpointer convert_thread_worker(gpointer p) {
	char destfilename[256], msg_bar[256];
	int i = 0;
	int indice;
	struct _convert_data *args = (struct _convert_data *) p;
	
	args->list = g_list_first (args->list);
	indice = args->start - 1;
	while(args->list) {
		char *sourcefilename;

		if (!get_thread_run()) {
			args->list = NULL;
			break;
		}
		sourcefilename = (char *)args->list->data;
		snprintf(destfilename, 255, "%s%05d.%s", destroot, ++indice, destsuf);
		++i;
		strncpy(sourcesuf, get_filename_ext(sourcefilename), SUFFIX_MAX_LEN);
		if (set_convflags_from_extension()) {
			char msg[512];
			siril_log_message("FILETYPE IS NOT SUPPORTED, CANNOT CONVERT: %s\n", sourcesuf);
			snprintf(msg, 511, "File extension '%s' is not supported.\n"
				"Verify that you typed the extension correctly.\n"
				"If so, you may need to install third-party software to enable "
				"this file type conversion, look at the README file.\n"
				"If the file type you are trying to load is listed in supported "
				"formats, you may notify the developpers that the extension you are "
				"trying to use should be recognized for this type.", sourcesuf);
			show_dialog(msg, "Warning", "gtk-dialog-warning");
		}
		else {
			if (convflags & CONVAVI){
				char newpath[256], newdest[256], newsource[256], *dirname;
				struct dirent *file;
				DIR *newdir;
				
				if (convert_film(sourcefilename)) {
					siril_log_message("Error during mplayer conversion. Aborted...\n");
					gdk_threads_add_idle(end_convert_idle, args);
					return NULL;
				}
				dirname = extract_path(sourcefilename);
				sprintf(newpath, "%s/%s", dirname, destroot);
				free(dirname);
				if((newdir = opendir(newpath)) == NULL){
					siril_log_message("Conversion: error opening working directory %s.\n", newpath);
					gdk_threads_add_idle(end_convert_idle, args);
					return NULL;
				}
				
				while ((file = readdir(newdir)) != NULL) {
					if (ends_with(file->d_name, destsuf)) continue;
					sprintf(newsource, "%s/%s", newpath, file->d_name);
					sprintf(newdest, "%s/seq_%s.%s", newpath, remove_ext_from_filename(file->d_name), destsuf);
					tofits(newsource, newdest);
					unlink(newsource);
				}
				closedir(newdir);
			}
			else
				tofits(sourcefilename, remove_ext_from_filename(destfilename));
		}
		char *name = strrchr(sourcefilename, '/');
		sprintf(msg_bar, "Converting %s...", name + 1);
		set_progress_bar_data(msg_bar, (double)i/((double)args->total));
		free(sourcefilename);
		args->list = g_list_next(args->list);
	}
	set_progress_bar_data(NULL, PROGRESS_DONE);

	gdk_threads_add_idle(end_convert_idle, args);
	return NULL;
}

// This idle function was not really required, but allows to properly join the thread.
static gboolean end_convert_idle(gpointer p) {
	struct _convert_data *args = (struct _convert_data *) p;
	struct timeval t_end;
	
	if (get_thread_run()) {
		// load the sequence of debayered images
		char *ppseqname = malloc(strlen(destroot) + 5);
		sprintf(ppseqname, "%s.seq", destroot);
		check_seq(0);
		update_sequences_list(ppseqname);
		free(ppseqname);
	}
	update_used_memory();
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_DONE);
	set_cursor_waiting(FALSE);
	gettimeofday(&t_end, NULL);
	show_time(args->t_start, t_end);
	stop_processing_thread();
	args->list = g_list_first(args->list);
	while (args->list) {
		g_free(args->list->data);
		args->list = g_list_next(args->list);
	}
	g_list_free (args->list);
	free(args);
	return FALSE;
}

int tofits(char *source, char *dest){
	char filename[256];
	int nbplan;
	fits *tmpfit = calloc(1, sizeof(fits));
	
/**********************************************************************
 * ***                     CONVERSION OF BMP                     **** *
 * *******************************************************************/

	if (convflags & CONVBMP){
		tmpfit->bitpix = USHORT_IMG;	// convert it to USHORT anyway
		nbplan = readbmp(source, tmpfit);
		switch (nbplan) {
			case -1:
				return 1;
			case 1:
				if(savefits(dest, tmpfit)){
					siril_log_message("tofits: savefit error, 1 plane\n");
				}
				break;
			case 3:
			case 4:
				if(convflags & CONV3X1){
					snprintf(filename, 255, "r_%s", dest);
					if (save1fits16(filename, tmpfit, RLAYER)) {
						siril_log_message("tofits: save1fit8 error, CONV3X1\n");
						return 1;
					}
					snprintf(filename, 255, "g_%s", dest);
					save1fits16(filename, tmpfit, GLAYER);
					snprintf(filename, 255, "b_%s", dest);
					save1fits16(filename, tmpfit, BLAYER);
				}
				else if(convflags & CONV1X3){
					//~ siril_log_message("tofits: CONV1X3\n");
					if(savefits(dest, tmpfit)){
						siril_log_message("tofits: savefit error, CONV1X3\n");
					}
				}
				else if(convflags & CONV1X1){
					if (save1fits16(dest, tmpfit, RLAYER)) {
						siril_log_message("tofits: save1fit8 error, CONV1X1\n");
					}
				}
				break;
			default:
				siril_log_message("Unrecognized BMP file %s, aborting (%d layers)\n",
						source, nbplan);
				return 1;
		}
	}
	
/**********************************************************************
 * ***                     CONVERSION OF PIC                     **** *
 * *******************************************************************/
	else if (convflags & CONVPIC){
		tmpfit->bitpix = USHORT_IMG;	// convert it to USHORT anyway but PIC are SHORT_IMG
		nbplan=readpic(source, tmpfit);	
		switch (nbplan){
			case -1:
				return 1;
			case 1:
				if(savefits(dest, tmpfit)){
					siril_log_message("tofits: savefit error, 1 plane\n");
				}
				break;
			case 3:
				if(convflags & CONV3X1){
					snprintf(filename, 255, "r_%s", dest);
					if (save1fits16(filename, tmpfit, RLAYER)) {
						siril_log_message("ERROR : Unable to convert the file %s.\n", filename);
						return 1;
					}
					snprintf(filename, 255, "g_%s", dest);
					save1fits16(filename, tmpfit, GLAYER);
					snprintf(filename, 255, "b_%s", dest);
					save1fits16(filename, tmpfit, BLAYER);;			
				}
				else if(convflags & CONV1X3){		
					savefits(dest, tmpfit);	
				}
				break;
			default:
				siril_log_message("Unrecognized PIC file %s, aborting (%d layers)\n",
						source, nbplan);
				return 1;
		}
	}
/**********************************************************************
 * ***                     CONVERSION OF TIFF                    **** *
 * *******************************************************************/	
#ifdef HAVE_LIBTIFF
	else if (convflags & CONVTIF){
		tmpfit->bitpix = USHORT_IMG;	// convert it to USHORT anyway
		nbplan=readtif(source, tmpfit);
		switch (nbplan){
			case -1:
				return 1;
			case 1:
				if(savefits(dest, tmpfit)){
					siril_log_message("tofits: savefit error, 1 plane\n");
				}
				break;	
			case 3:		
				if(convflags & CONV3X1){
					snprintf(filename, 255, "r_%s", dest);
					if (save1fits16(filename, tmpfit, RLAYER)) {
						siril_log_message("ERROR : Unable to convert the file %s\n", filename);
						return 1;
					}
					snprintf(filename, 255, "g_%s", dest);
					save1fits16(filename, tmpfit, GLAYER);
					snprintf(filename, 255, "b_%s", dest);
					save1fits16(filename, tmpfit, BLAYER);			
				}
				else if(convflags & CONV1X3){		
					savefits(dest, tmpfit);	
				}
				break;
			default:
				siril_log_message("Unrecognized TIFF file %s, aborting (%d layers)\n",
						source, nbplan);
				return 1;			
			}
	}
#endif

/**********************************************************************
 * ***                     CONVERSION OF JPG                     **** *
 * *******************************************************************/
#ifdef HAVE_LIBJPEG
	else if (convflags & CONVJPG){
		tmpfit->bitpix = USHORT_IMG;	// convert it to USHORT anyway
		nbplan=readjpg(source, tmpfit);	
		switch (nbplan){
			case -1:
				return 1;
			case 1:
				if(savefits(dest, tmpfit)){
					siril_log_message("tofits: savefit error, 1 plane\n");
				}
				break;	
			case 3:	
				if(convflags & CONV3X1){
					snprintf(filename, 255, "r_%s", dest);
					if (save1fits16(filename, tmpfit, RLAYER)) {
						siril_log_message("ERROR : Unable to convert the file %s\n", filename);
						return 1;
					}
					snprintf(filename, 255, "g_%s", dest);
					save1fits16(filename, tmpfit, GLAYER);
					snprintf(filename, 255, "b_%s", dest);
					save1fits16(filename, tmpfit, BLAYER);;			
				}
				else if(convflags & CONV1X3){		
					savefits(dest, tmpfit);	
				}
				break;
			default:
				siril_log_message("Unrecognized JPEG file %s, aborting (%d layers)\n",
						source, nbplan);
				return 1;			
			}		
	}
#endif

/**********************************************************************
 * ***                     CONVERSION OF PNG                     **** *
 * *******************************************************************/
#ifdef HAVE_LIBPNG
	else if (convflags & CONVPNG){
		tmpfit->bitpix = USHORT_IMG;	// convert it to USHORT anyway
		nbplan=readpng(source, tmpfit);
		switch(nbplan){	
			case -1:
				return 1;
			case 1:
				if(savefits(dest, tmpfit)){
					siril_log_message("tofits: savefit error, 1 plane\n");
				}
				break;	
			case 3:	
				if(convflags & CONV3X1){
					snprintf(filename, 255, "r_%s", dest);
					if (save1fits16(filename, tmpfit, RLAYER)) {
						siril_log_message("ERROR : Unable to convert the file %s\n", filename);
						return 1;
					}
					snprintf(filename, 255, "g_%s", dest);
					save1fits16(filename, tmpfit, GLAYER);
					snprintf(filename, 255, "b_%s", dest);
					save1fits16(filename, tmpfit, BLAYER);;			
				}
				else if(convflags & CONV1X3){		
					savefits(dest, tmpfit);	
				}
				break;
			default:
				siril_log_message("Unrecognized PNG file %s, aborting (%d layers)\n",
						source, nbplan);
				return 1;			
			}					
	}
#endif

/**********************************************************************
 * ***                     CONVERSION OF RAW                     **** *
 * *******************************************************************/
#ifdef HAVE_LIBRAW
	else if (convflags & CONVRAW) {
		open_raw_files(source, tmpfit, com.raw_set.cfa);
		if(convflags & CONV3X1){
			if (com.raw_set.cfa) {
				siril_log_message("Cannot convert the B&W file into 3 channels\n");
				clearfits(tmpfit);
				return 1;
			}
			snprintf(filename, 255, "r_%s", dest);
			if (save1fits16(filename, tmpfit, RLAYER)) {
				siril_log_message("ERROR : Unable to convert the file %s\n", filename);
				return 1;
			}
			snprintf(filename, 255, "g_%s", dest);
			save1fits16(filename,tmpfit,GLAYER);
			snprintf(filename, 255, "b_%s", dest);
			save1fits16(filename, tmpfit, BLAYER);;			
		}
		else if(convflags & CONV1X3){		
			savefits(dest, tmpfit);
		}
	}
#endif

/**********************************************************************
 * ***                   CONVERSION OF NetPBM                    **** *
 * *******************************************************************/
	else if (convflags & CONVPNM){
		tmpfit->bitpix = USHORT_IMG;	// convert it to USHORT anyway
		nbplan = import_pnm_to_fits(source, tmpfit);
		switch (nbplan) {
			case -1:
				return 1;
			case 1:
				if(savefits(dest, tmpfit)){
					siril_log_message("tofits: savefit error, 1 plane\n");
				}
				break;
			case 3:
				if(convflags & CONV3X1){
					snprintf(filename, 255, "r_%s", dest);
					if (save1fits16(filename, tmpfit, RLAYER)) {
						siril_log_message("tofits: save1fit8 error, CONV3X1\n");
						return 1;
					}
					snprintf(filename ,255, "g_%s", dest);
					save1fits16(filename, tmpfit, GLAYER);
					snprintf(filename, 255, "b_%s", dest);
					save1fits16(filename, tmpfit, BLAYER);
				}
				else if(convflags & CONV1X3){
					//~ siril_log_message("tofits: CONV1X3\n");
					if(savefits(dest, tmpfit)){
						siril_log_message("tofits: savefit error, CONV1X3\n");
					}
				}
				else if(convflags & CONV1X1){
					if (save1fits16(dest, tmpfit, RLAYER)) {
						siril_log_message("tofits: save1fit8 error, CONV1X1\n");
					}
				}
				break;
			default:
				siril_log_message("Unrecognized NetPNM file %s, aborting (%d layers)\n",
						source, nbplan);
				return 1;
		}
	}

/**********************************************************************
 * ***                   CONVERSION OF FITS CFA                  **** *
 * *******************************************************************/
/* Siril's FITS are stored bottom to top, debayering will throw 
 * wrong results. So before demosacaing we need to transforme the image
 * with fits_flip_top_to_bottom() function */
	else if (convflags & CONVCFA) {
		readfits(source, tmpfit, NULL);
		fits_flip_top_to_bottom(tmpfit);
		siril_log_message("Filter Pattern: %s\n", filter_pattern[com.raw_set.bayer_pattern]);
		if (!debayer(tmpfit, com.raw_set.bayer_inter)) {
			fits_flip_top_to_bottom(tmpfit);
			savefits(dest, tmpfit);
		} else {
			siril_log_message("Cannot perform debayering\n");
		}
	}

/**********************************************************************
 * ***                     CONVERSION OF FITS                    **** *
 * *******************************************************************/
	else if (convflags & CONVFIT) {
		readfits(source, tmpfit, NULL);
		if (tmpfit->naxes[2] == 3 && (convflags & CONV3X1)){
			snprintf(filename, 255, "r_%s", dest);
			if (save1fits16(filename, tmpfit, RLAYER)) {
				return 1;
			}
			snprintf(filename ,255, "g_%s", dest);
			save1fits16(filename, tmpfit, GLAYER);
			snprintf(filename, 255, "b_%s", dest);
			save1fits16(filename, tmpfit, BLAYER);
		}
		else {
			savefits(dest, tmpfit);
		}
	}

	clearfits(tmpfit);
	return 0;
}
