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

#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include <dirent.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_histogram.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "siril.h"
#include "command.h"
#include "proto.h"
#include "conversion.h"
#include "callbacks.h"
#include "colors.h"
#include "PSF.h"
#include "PSF_list.h"
#include "star_finder.h"
#include "Def_Math.h"
#include "Def_Wavelet.h"
#include "histogram.h"
#include "single_image.h"
#include "gradient.h"
#include "fft.h"

#ifdef HAVE_OPENCV
#include "opencv.h"
#endif

static char *word[MAX_COMMAND_WORDS];

command commande[] = {
	/* name,	nbarg,	usage,			function pointer */
	{"addmax",	1,	"addmax filename",	process_addmax},
	
	{"bg", 0, "bg", process_bg},
	{"bgnoise", 0, "bgnoise", process_bgnoise},
	
	{"cd", 1, "cd directory (define the working directory)", process_cd},
	{"clearstar", 0, "clearstar", process_clearstar},
	{"contrast", 0, "contrast", process_contrast},
	{"crop", 0, "crop [x y width height]", process_crop}, 
//	{"crop2", 3, "crop2 genname outname number [x y width height]", process_crop2}, 

	{"ddp", 3, "ddp level coef sigma", process_ddp}, 
	
	{"entropy", 0, "entropy", process_entropy},
	{"exit", 0, "exit", process_exit},
	{"extract", 1, "extract NbPlans", process_extract},
	
	{"fdiv", 2, "fdiv filename scalar", process_fdiv},
	{"fill", 1, "fill value", process_fill},
	{"fill2", 1, "fill2 value [x y width height]", process_fill2},
	{"findstar", 0, "findstar", process_findstar},
	//~ {"find_hot", 2, "find_hot file threshold", process_findhot},
	{"fmedian", 2, "fmedian ksize modulation", process_fmedian},
	{"fftd", 2, "fftd magnitude phase", process_fft},
	{"ffti", 2, "ffti magnitude phase", process_fft},
	{"fixbanding", 2, "fixbanding amount", process_fixbanding},
	
	{"gauss", 1, "gauss sigma ", process_gauss},	
	//~ {"gauss2", 1, "gauss sigma", process_gauss2},

	{"help", 0, "help", process_help},	
	{"histo", 1, "histo layer (layer=0, 1, 2 with 0: red, 1: green, 2: blue)", process_histo},
	
	/* i*** commands oper filename and curent image */
	{"iadd", 1, "add filename", process_imoper}, 
	{"idiv", 1, "idiv filename", process_imoper},
	{"imul", 1, "imul filename", process_imoper}, 
	{"isub", 1, "isub filename", process_imoper},
	
	{"load", 1, "load filename.[ext]", process_load}, 
	// specific loads are not required, but could be used to force the
	// extension to a higher priority in case two files with same basename
	// exist (stat_file() manages that priority order for now).
	{"log", 0, "log ", process_log}, /* logarifies current image */
	{"ls", 0, "ls ", process_ls},
	
	{"mirrorx", 0, "mirrorx", process_mirrorx},
	{"mirrory", 0, "mirrory", process_mirrory},
	
	{"new", 2, "new width height nb_layers", process_new},
	{"nozero", 1, "nozero level (replaces null values by level)", process_nozero}, /* replaces null values by level */
	
	{"offset", 1, "offset value", process_offset},
	
	{"psf", 0, "psf", process_psf},
	
#ifdef HAVE_OPENCV
	{"resample", 1, "resample factor", process_resample},
#endif	
	{"rmgreen", 1, "rmgreen type", process_scnr},
#ifdef HAVE_OPENCV
	{"rotate", 1, "rotate angle", process_rotate},
#endif
	{"rotatePi", 0, "rotatePi", process_rotatepi},
	
	{"satu", 1, "satu coeff ", process_satu}, 
	{"save", 1, "save filename (save current image in fit)", process_save}, 
	{"savebmp", 1, "savebmp filename (save display image in bmp)", process_savebmp}, 
#ifdef HAVE_LIBJPEG
	{"savejpg", 1, "savejpg filename [quality] (save current display in jpg)", process_savejpg},
#endif
	{"savepnm", 1, "savepnm filename (save current image in Netpbm)", process_savepnm},
#ifdef HAVE_LIBTIFF
	{"savetif", 1, "savetif filename (save current image in tif 16bits)", process_savetif},
	{"savetif8", 1, "savetif8 filename (save current image in tif 8bits)", process_savetif},
#endif
	{"split", 3, "split R G B", process_split},
	{"stat", 0, "stat", process_stat},
	
#ifdef _OPENMP
	{"setcpu", 1, "setcpu number", process_set_cpu},
#endif
	{"threshlo", 1, "threshlo level", process_threshlo},
	{"threshhi", 1, "threshi level", process_threshhi}, 
	{"thresh", 2, "thresh hi lo (threshes hi and lo)", process_thresh}, /* threshes hi and lo */
	
	/* unsharp masking of current image or genname sequence */
	{"unsharp", 2, "unsharp sigma multi", process_unsharp},
//	{"unsharp2", 5, "unsharp2 sigma multi src dest number", process_unsharp2},

	{"visu", 2, "visu low high", process_visu},
	
	/* wavelet transform in nbr_plan plans */ 
	{"wavelet", 1, "wavelet nbr_plan type (1=linear 2=spline)", process_wavelet},
	/* reconstruct from wavelet transform and weighs plans with c1, c2, c3... */ 
	{"wrecons", 2, "wrecons c1 c2 c3 ...", process_wrecons},
	
	{"",0,"",0}
};

int process_load(int nb){
	char filename[256];
	int retval, i;
	
	strncpy(filename, word[1], 250);
	filename[250] = '\0';	
	
	for (i=1; i<nb-1; ++i){
		strcat(filename, " ");
		strcat(filename, word[i+1]);
	}
	expand_home_in_filename(filename, 256);
	retval = open_single_image(filename);
	return retval;
}

int process_satu(int nb){
	if (get_thread_run()) {
		siril_log_message("Another task is already in progress, ignoring new request.\n");
		return 1;
	}
	struct enhance_saturation_data *args = malloc(sizeof(struct enhance_saturation_data));
	
	args->coeff = atof(word[1]);
	if (args->coeff == 0.0) args->coeff = 1.0;

	args->fit = &gfit;
	args->h_min = 0.0;
	args->h_max = 360.0;
	args->preserve = TRUE;
	set_cursor_waiting(TRUE);
	start_in_new_thread(enhance_saturation, args);
	
	return 0;
}

int process_save(int nb){
	char filename[256];
	
	if (sequence_is_loaded() && !single_image_is_loaded()) {
		gfit.hi = com.seq.layers[RLAYER].hi;
		gfit.lo = com.seq.layers[RLAYER].lo;
	}
	else {
		gfit.hi = com.uniq->layers[RLAYER].hi;
		gfit.lo = com.uniq->layers[RLAYER].lo;
	}


	sprintf(filename, "%s", word[1]);
	set_cursor_waiting(TRUE);
	savefits(filename, &(gfit));
	set_cursor_waiting(FALSE);
	return 0;
}

int process_savebmp(int nb){
	char filename[256];
	
	sprintf(filename, "%s", strcat(word[1], ".bmp"));
	set_cursor_waiting(TRUE);
	savebmp(filename, &(gfit));
	set_cursor_waiting(FALSE);
	return 0;
}

#ifdef HAVE_LIBJPEG
int process_savejpg(int nb){
	char filename[256];
	int quality = 100;
	
	if ((nb == 3) && atoi(word[2]) <= 100 && atoi(word[2]) > 0)
		quality=atoi(word[2]);
	strcpy(filename, word[1]);
	strcat(filename, ".jpg");
	set_cursor_waiting(TRUE);
	savejpg(filename, &gfit, quality);
	set_cursor_waiting(FALSE);
	return 0;
}
#endif

#ifdef HAVE_LIBTIFF
int process_savetif(int nb){
	char filename[256];
	uint16 bitspersample = 16;
	
	if (strcasecmp(word[0],"savetif8")==0) bitspersample=8;
	sprintf(filename,"%s", strcat(word[1],".tif"));
	set_cursor_waiting(TRUE);
	savetif(filename, &gfit, bitspersample);
	set_cursor_waiting(FALSE);
	return 0;
}
#endif

int process_savepnm(int nb){
	char filename[256];
	
	switch(gfit.naxes[2]){
		case 1:
			sprintf(filename,"%s", strcat(word[1],".pgm"));
				set_cursor_waiting(TRUE);
				savepgm(filename, &(gfit));
				set_cursor_waiting(FALSE);
		break;
		case 3:
			sprintf(filename,"%s", strcat(word[1],".ppm"));
			set_cursor_waiting(TRUE);
			saveppm(filename, &(gfit));
			set_cursor_waiting(FALSE);
		break;
		
		default:	/*Should not happend */
			return 1;
		}
	return 0;	
}

int process_imoper(int nb){
	clearfits(&(wfit[4]));
	readfits(word[1], &(wfit[4]), NULL);
	imoper(&gfit, &wfit[4], word[0][1]);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_addmax(int nb){
	clearfits(&(wfit[4]));
	readfits(word[1], &(wfit[4]), NULL);
	if (addmax(&gfit, &wfit[4])==0) {
		adjust_cutoff_from_updated_gfit();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
	}
	return 0;
}

int process_fdiv(int nb){
	// combines an image division and a scalar multiplication.
	float norm;

	norm = atof(word[2]);
	clearfits(&(wfit[4]));
	readfits(word[1], &(wfit[4]), NULL);
	fdiv(&gfit,&wfit[4], norm);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_entropy(int nb){
	rectangle area;
	if (com.selection.w > 0 && com.selection.h > 0)
		memcpy(&area, &com.selection, sizeof(rectangle));
	else {
		area.x = area.y = 0;
		area.h = gfit.ry;
		area.w = gfit.rx;
	}
	double ent = entropy(&gfit, com.cvport, &area, NULL);
	siril_log_message("entropy: %10.2f\n", ent);
	return 0;
}

int process_gauss(int nb){
	unsharp(&(gfit), atof(word[1]), (double)0, TRUE);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_gauss2(int nb){
	int i;
	
	for (i=0; i<com.seq.number; ++i) {
		printf("test\n");
	}
	return 0;
}

int process_unsharp(int nb){
	unsharp(&(gfit), atof(word[1]), atof(word[2]), TRUE);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_crop(int nb){
	rectangle area;
	if (!com.drawn || com.drawing){		// TODO: what's that test?
		if (nb==5){
			if (atoi(word[1])<0 || atoi(word[2])<0){
				siril_log_message("Crop: x and y must be positive values.\n");
				return 1;
			}			
			if (atoi(word[3])<=0 || atoi(word[4])<=0){
				siril_log_message("Crop: width and height must be greater than 0.\n");
				return 1;
			}
			if (atoi(word[3])>gfit.rx || atoi(word[4])>gfit.ry){
				siril_log_message("Crop: width and height, respectively, must be less than %d and %d.\n", gfit.rx,gfit.ry);
				return 1;
			}
			area.x = atoi(word[1]);
			area.y = atoi(word[2]);
			area.w = atoi(word[3]);
			area.h = atoi(word[4]);
		}
		else {
			siril_log_message("Crop: select a region or provide x,y,width,height\n");
			return 1;
		}
	} else {
		memcpy(&area, &com.selection, sizeof(rectangle));
	}

	crop(&gfit, &area);
	delete_selected_area();
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	update_used_memory();
	return 0;
}

int process_cd(int nb){
	char filename[256];
	int retval, i;
	
	strncpy(filename, word[1], 250);
	filename[250] = '\0';
	
	for (i=1; i<nb-1; ++i){
		strcat(filename, " ");
		strcat(filename, word[i+1]);
	}
	
	expand_home_in_filename(filename, 256);
	if (!(retval=changedir(filename)))
		writeinitfile();
	return retval;
}

int process_wrecons(int nb){
	int i;
	float coef[7];
	char *File_Name_Transform[3] = {"r_rawdata.wave", "g_rawdata.wave", "b_rawdata.wave"}, *dir[3];
	const char *tmpdir;
	int nb_chan = gfit.naxes[2];
	
	assert(nb_chan == 1 || nb_chan == 3);
		
 	tmpdir = g_get_tmp_dir();
 	
	for (i=0;i<nb-1;++i){
		coef[i] = atof(word[i+1]);
	}

	for (i=0; i < nb_chan; i++) {
		dir[i] = malloc(strlen(tmpdir) + strlen(File_Name_Transform[i]) + 2);
		strcpy(dir[i], tmpdir);
		strcat(dir[i], "/");
		strcat(dir[i], File_Name_Transform[i]);
		wavelet_reconstruct_file (dir[i], coef, gfit.pdata[i]);
		free(dir[i]);
	}

	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}	

int process_wavelet(int nb){
	char *File_Name_Transform[3] = {"r_rawdata.wave", "g_rawdata.wave", "b_rawdata.wave"}, *dir[3];
	const char* tmpdir;
	int Type_Transform, Nbr_Plan, maxplan, mins, chan, nb_chan;
	float *Imag;
	
 	tmpdir = g_get_tmp_dir();
	
	Nbr_Plan = atoi(word[1]);
	Type_Transform = atoi(word[2]);
	
	nb_chan = gfit.naxes[2];
	assert(nb_chan <= 3);

	mins = min (gfit.rx, gfit.ry);
	maxplan = log(mins) / log(2) - 2;

	if ( Nbr_Plan > maxplan ){
		siril_log_message("Wavelet: maximum number of plans for this image size is %d\n", 
				maxplan);
		return 1;
	}

	if(Type_Transform != TO_PAVE_LINEAR && Type_Transform !=TO_PAVE_BSPLINE){
		siril_log_message("Wavelet: type must be %d or %d\n",TO_PAVE_LINEAR,TO_PAVE_BSPLINE);
		return 1;
	}

	Imag = f_vector_alloc (gfit.rx * gfit.ry);
	
	for (chan = 0; chan < nb_chan; chan++) {
		dir[chan] = malloc(strlen(tmpdir) + strlen(File_Name_Transform[chan]) + 2);
		strcpy(dir[chan], tmpdir);
		strcat(dir[chan], "/");
		strcat(dir[chan], File_Name_Transform[chan]);
		wavelet_transform_file (Imag, dir[chan], Type_Transform, Nbr_Plan, gfit.pdata[chan]);
		free(dir[chan]);
	}
	
	free (Imag);
	return 0;
}

int process_log(int nb){
	loglut(&gfit, LOG);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_ls(int nb){
	DIR* ptdir;
	struct dirent* entry;
	char filename[256];
	
	filename[0]='\0';
	/* If a path is given in argument */
	if (nb>1){
		if (word[1][0]!='\0'){	
			/* Absolute path */
			if(word[1][0]=='/' || word[1][0]=='~'){
				strncpy(filename, word[1], 250);
				filename[250] = '\0';
				expand_home_in_filename(filename, 256);
			}
			/* Relative path */
			else {
				strcpy(filename, com.wd);
				strcat(filename, "/");
				strcat(filename, word[1]);
			}	
			ptdir = opendir(filename);
		}
		/* Should not happend */
		else {
			printf("Cannot list files in %s\n", filename);
			return 1;
		}
	}
	/* No paths are given in argument */
	else {
		if (!com.wd) {
			siril_log_message("Cannot list files, set working directory first.\n");
			return 1;
		}
		ptdir = opendir(com.wd);
	}
	if (!ptdir) {
		siril_log_message("Siril cannot open the directory.\n");
		return 1;
	}
	/* List the entries */
	while ((entry = readdir(ptdir)) != NULL) {
		struct stat entrystat;
		char file_path[256];
		if (entry->d_name[0] == '.')
			continue;	/* no hidden files */
		if (filename[0] != '\0')
			sprintf(file_path, "%s/%s", filename, entry->d_name);
		else
			sprintf(file_path, "%s", entry->d_name);
		if (lstat(file_path, &entrystat)) {
			perror("stat");
			break;		
		}
		if (S_ISLNK(entrystat.st_mode)) {
			siril_log_color_message("Link: %s\n", "bold", entry->d_name);
			continue;
		}
		if (S_ISDIR(entrystat.st_mode)) {
			siril_log_color_message("Directory: %s\n", "green", entry->d_name);
			continue;
		}
		int fnlen = strlen(entry->d_name);
		if (fnlen < 5) continue;
		int extlen = strlen(get_filename_ext(entry->d_name));
		if (get_type_for_extension_name(entry->d_name + fnlen - (extlen + 1)) != TYPEUNDEF){
			if ((get_type_for_extension_name(entry->d_name + fnlen - (extlen + 1)) == TYPEAVI) ||
					(get_type_for_extension_name(entry->d_name + fnlen - 4) == TYPESER))
				siril_log_color_message("Video: %s\n", "salmon", entry->d_name);
			else
				siril_log_color_message("Image: %s\n", "red", entry->d_name);
		}
		// RAW files are not listed with the above filter
		if (!strncasecmp(entry->d_name + fnlen - 4, ".seq", 4))
			siril_log_color_message("Sequence: %s\n", "blue", entry->d_name);
	}
	siril_log_message("********* END OF THE LIST *********\n");
	closedir(ptdir);
	return 0;
}

int	process_mirrorx(int nb){
	mirrorx(&gfit, TRUE);
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int	process_mirrory(int nb){
	mirrory(&gfit, TRUE);
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

#ifdef HAVE_OPENCV
int process_resample(int nb) {
	double factor = atof(word[1]);
	if (factor > 5.0) {
		siril_log_message("The scaling factor must be less than 5.0\n");
		return 1;
	}
	int toX = round_to_int(factor * gfit.rx);
	int toY = round_to_int(factor * gfit.ry);
	
	set_cursor_waiting(TRUE);
	verbose_resize_gaussian(&gfit, toX, toY, 1);
	update_used_memory();
	adjust_vport_size_to_image();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();	
	set_cursor_waiting(FALSE);
	return 0;
}


int process_rotate(int nb) {
	double degree;
	
	set_cursor_waiting(TRUE);
	degree = atof(word[1]);
	verbose_rotate_image(&gfit, degree, 1, 1);	//INTER_LINEAR
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	return 0;
}
#endif

int process_rotatepi(int nb){
#ifdef HAVE_OPENCV
	verbose_rotate_image(&gfit, 180.0, 1, 1);
#else
	fits_rotate_pi(&gfit);
#endif
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;	
}

int process_psf(int nb){
	int layer = match_drawing_area_widget(com.vport[com.cvport], FALSE);
	if (layer != -1) {
		if (!com.drawn || com.drawing)
			return 1;
		if (com.selection.w > 300 || com.selection.h > 300){
			siril_log_message("Current selection is too large. To determine the PSF, please make a selection around a star.\n");
			return 1;
		}
		fitted_PSF *result = get_Minimisation(&gfit, layer, &com.selection);
		if (result) {
			DisplayResult(result, &com.selection);
			free(result);
		}
	}
	return 0;
}

int process_bg(int nb){
	WORD bg = round_to_WORD(background(&gfit, -1, &com.selection));
	siril_log_message("Background value: %d\n", bg);
	return 0;
}

int process_bgnoise(int nb){
	if (get_thread_run()) {
		siril_log_message(
				"Another task is already in progress, ignoring new request.\n");
		return 1;
	}

	struct noise_data *args = malloc(sizeof(struct noise_data));

	set_cursor_waiting(TRUE);
	control_window_switch_to_tab(OUTPUT_LOGS);

	args->fit = &gfit;
	args->verbose = TRUE;
	memset(args->bgnoise, 0.0, sizeof(double[3]));
	start_in_new_thread(noise, args);
	return 0;
}

int process_histo(int nb){
	size_t i;
	int nlayer = atoi(word[1]);
	char* clayer;
	char name [20];
	
	if (nlayer>3 || nlayer <0)
		return 1;
	gsl_histogram* histo = computeHisto(&gfit, nlayer);
	if (!isrgb(&gfit)) clayer = strdup("bw");		//if B&W
	else clayer = vport_number_to_name(nlayer);
	snprintf(name, 20, "histo_%s.dat",clayer);

	FILE *f = fopen(name, "w");

	if (f == NULL) {
		free(clayer);
		return 1;
	}
	for (i=0; i < USHRT_MAX + 1; i++)
		fprintf(f, "%d %d\n",(unsigned int)i, (unsigned int)gsl_histogram_get (histo, i));
	fclose(f);
	gsl_histogram_free(histo);
	siril_log_message("The file %s has been created for the %s layer.\n", name, clayer);
	return 0;
}

int process_thresh(int nb){
	int lo, hi;

	lo = atoi(word[1]);
	hi = atoi(word[2]);
	threshlo(&gfit, lo);
	threshhi(&gfit, hi);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_threshlo(int nb){
	int lo;

	lo = atoi(word[1]);
	threshlo(&gfit, lo);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_threshhi(int nb){
	int hi;

	hi = atoi(word[1]);
	threshhi(&gfit, hi);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_nozero(int nb){
	int level;

	level = atoi(word[1]);
	nozero(&gfit, level);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_ddp(int nb){
	// combines an image division and a scalar multiplication.
	float coeff, sigma;
	unsigned level;

	level = atoi(word[1]);
	coeff = atof(word[2]);
	sigma = atof(word[3]);
	ddp(&gfit, level, coeff, sigma);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_new(int nb){
	int width, height, layers;
	
	width = atof(word[1]);
	height = atof(word[2]);
	layers = atoi(word[3]);
	if (layers != 1 && layers != 3) {
		siril_log_message("Number of layers MUST be 1 or 3\n");
		return 1;
	}
	if (!height || !width) return 1;

	close_single_image();

	new_fit_image(&gfit, width, height, layers);

	open_single_image_from_gfit(strdup("new empty image"));
	return 0;
}

int process_visu(int nb){
	int low, high;
	
	low = atoi(word[1]);
	high = atoi(word[2]);
	if ((high>USHRT_MAX) || (low<0)){
		siril_log_message("Values must be positive and less than %d.\n", USHRT_MAX);
		return 1;		
	}
	visu(&gfit, low, high);
	return 0;
}

int process_fill2(int nb){
	int level=atoi(word[1]);
	rectangle area;
	if (!com.drawn || com.drawing){		// TODO: what's that test?
		if (nb==6){
			area.x = atoi(word[2]);
			area.y = atoi(word[3]);
			area.w = atoi(word[4]);
			area.h = atoi(word[5]);
		}
		else {
			siril_log_message("Fill2: select a region or provide x,y,width,height\n");
			return 1;
		}
	} else {
		memcpy(&area, &com.selection, sizeof(rectangle));
	}
	fill(&gfit, level, &area);
	area.x = gfit.rx - area.x - area.w;
	area.y = gfit.ry - area.y - area.h;
	fill(&gfit, level, &area);
	redraw(com.cvport, REMAP_ALL);
	return 0;
}

int process_findstar(int nb){
	int layer = RLAYER;
	starFinder sf;

	memset(&sf, 0, sizeof(starFinder));
	
	if (!single_image_is_loaded()) return 0;
	if (isrgb(&gfit)) layer = GLAYER;
	delete_selected_area();
	com.stars = peaker(&gfit, layer, &sf);
	refresh_stars_list(com.stars);
	return 0;
}

/* FIXME */
int process_findhot(int nb){
	double sigma = atof(word[2]);
	int count;
	
	count = find_hot_pixels(&gfit, sigma, word[1]);
	siril_log_message("Number of hot pixels: %d\n", count);
	//~ adjust_cutoff_from_updated_gfit();
	//~ redraw(com.cvport, REMAP_ALL);
	//~ redraw_previews();
	return 0;
}

int process_fmedian(int nb){
	if (get_thread_run()) {
		siril_log_message("Another task is already in progress, ignoring new request.\n");
		return 1;
	}
	
	struct median_filter_data *args = malloc(sizeof(struct median_filter_data));	
	args->ksize = atoi(word[1]);
	args->amount = atof(word[2]);
	args->iterations = 1;
	
	if (!(args->ksize & 1) || args->ksize < 2) {
		siril_log_message("The size of the kernel MUST be odd and greater than 1.\n");
		free(args);
		return 1;
	}
	if (args->amount < 0.0 || args->amount > 1.0) {
		siril_log_message("Modulation value MUST be between 0 and 1\n");
		free(args);
		return 1;
	}
	args->fit = &gfit;

	set_cursor_waiting(TRUE);
	start_in_new_thread(median_filter, args);
	
	return 0;
}

int process_clearstar(int nb){
	clear_stars_list();
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_NONE);
	redraw_previews();
	return 0;
}

int process_contrast(int nb){
	int layer;
	double result[gfit.naxes[2]], value=0;
	
	for (layer = 0; layer < gfit.naxes[2]; layer++)
		result[layer] = contrast(&gfit, layer);
	for (layer = 0; layer < gfit.naxes[2]; layer++)
		value += result[layer];
	value /= gfit.naxes[2];
	
	siril_log_message("Contrast: %lf\n", value);
	return 0;
}

int process_fill(int nb){	
	int level;
	rectangle area;
	
	if (!com.drawn || com.drawing){		// TODO: what's that test?
		if (nb==6){
			area.x = atoi(word[2]);
			area.y = atoi(word[3]);
			area.w = atoi(word[4]);
			area.h = atoi(word[5]);
		}
		else {
			area.w = gfit.rx; area.h = gfit.ry;
			area.x = 0; area.y = 0;
		}
	} else {
		memcpy(&area, &com.selection, sizeof(rectangle));
	}
	level = atoi(word[1]);
	fill(&gfit, level, &area);
	redraw(com.cvport, REMAP_ALL);
	return 0;
}

int process_offset(int nb){
	int level;
	
	level = atoi(word[1]);
	off(&gfit, level);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

/* The version in command line is a minimal version
 * Only neutral type are available (no amount needed), 
 * then we always preserve the lightness */
int process_scnr(int nb){
	if (get_thread_run()) {
		siril_log_message("Another task is already in progress, ignoring new request.\n");
		return 1;
	}
	struct scnr_data *args = malloc(sizeof(struct scnr_data));
	
	args->type = atoi(word[1]);
	args->fit = &gfit;
	args->amount = 0.0;
	args->preserve = TRUE;
	set_cursor_waiting(TRUE);
	start_in_new_thread(scnr, args);

	return 0;
}

int process_fft(int nb){
	if (sequence_is_loaded()) {
		siril_log_message("FFT does not work with sequences\n");
		return 1;
	}
	if (get_thread_run()) {
		siril_log_message("Another task is already in progress, ignoring new request.\n");
		return 1;
	}
	struct fft_data *args = malloc(sizeof(struct fft_data));
	
	args->fit = &gfit;
	args->type = strdup(word[0]);
	args->modulus = strdup(word[1]);
	args->phase = strdup(word[2]);
	args->type_order = 0;
	
	set_cursor_waiting(TRUE);
	start_in_new_thread(fourier_transform, args);
	
	return 0;
}

int process_fixbanding(int nb) {
	if (get_thread_run()) {
		siril_log_message("Another task is already in progress, ignoring new request.\n");
		return 1;
	}
	
	struct banding_data *args = malloc(sizeof(struct banding_data));

	args->amount = atof(word[1]);
	args->sigma = atof(word[2]);
	args->protect_highlights = TRUE;
	args->fit = &gfit;

	set_cursor_waiting(TRUE);
	start_in_new_thread(BandingEngine, args);
	
	return 0;
}

int process_split(int nb){
	char R[256], G[256], B[256];
	
	if (!isrgb(&gfit)) {
		siril_log_message("Siril cannot split layers. Make sure your image is in RGB mode.\n");
		return 1;
	}
	sprintf(R, "%s.fit", word[1]);
	sprintf(G, "%s.fit", word[2]);
	sprintf(B, "%s.fit", word[3]);
	save1fits16(R, &gfit, RLAYER);
	save1fits16(G, &gfit, GLAYER);
	save1fits16(B, &gfit, BLAYER);
	return 0;
}

int process_stat(int nb){
	int nplane = gfit.naxes[2];
	int layer;
	
	for (layer=0;layer<nplane;layer++){
		imstats* stat = statistics(&gfit, layer, &com.selection);
		siril_log_message("%s layer: Mean: %0.1lf, Median: %0.1lf, Sigma: %0.1lf, Min: %0.1lf, Max: %0.1lf\n",
					stat->layername, stat->mean, stat->median, stat->sigma, stat->min, stat->max);
		free(stat);
		stat = NULL;
	}
	return 0;
}

#ifdef _OPENMP
int process_set_cpu(int nb){
	int proc_in, proc_out, proc_max;

	proc_in = atoi(word[1]);
	proc_max = omp_get_num_procs();
	if (proc_in > proc_max || proc_in < 1) {
		siril_log_message(
				"Number of logical processor MUST be greater than 0 and lower or equal to %d.\n", proc_max);
		return 1;
	}
	omp_set_num_threads(proc_in);

#pragma omp parallel
	{
		proc_out = omp_get_num_threads();
	}
	siril_log_message("Using now %d logical processor\n", proc_out);
	com.max_thread = proc_out;
	update_spinCPU(0);

	return 0;
}
#endif

int process_help(int nb){
	command *current = commande;
	siril_log_message("********* LIST OF AVAILABLE COMMANDS *********\n");
	while(current->process){
		siril_log_message("%s\n", current->usage);
		current++;
	}
	siril_log_message("********* END OF THE LIST *********\n");
	return 0;
}

int process_exit(int nb){
	exit(0);
}

int process_extract(int nb) {
	int Nbr_Plan, maxplan, mins, i;
	
	Nbr_Plan = atoi(word[1]);

	mins = min (gfit.rx, gfit.ry);
	maxplan = log(mins) / log(2) - 2;

	if ( Nbr_Plan > maxplan ){
		siril_log_message("Wavelet: maximum number of plans for this image size is %d\n", 
				maxplan);
		return 1;
	}
	fits *fit = calloc(1, sizeof(fits));
	copyfits(&gfit, fit, CP_ALLOC | CP_COPYA | CP_FORMAT, 0);
	
	for (i=0; i < Nbr_Plan; i++) {
		char filename[256];
		
		sprintf(filename, "layer%02d", i);
		get_wavelet_layers(fit, Nbr_Plan, i, TO_PAVE_BSPLINE, -1);
		savefits(filename, fit);
	}
	clearfits(fit);
	update_used_memory();
	return 0;
}

int processcommand(const char *line) {
	int i = 0, wordnb = 0, len = strlen(line);
	char *myline;
	char string_starter = '\0';	// quotes don't split words on spaces
	word[0] = NULL;
	if (line[0] == '\0' || line[0] == '\n') return 0;
	myline = strdup(line);
	do {
		while (i<len && isblank(myline[i])) i++;
		if (myline[i] == '"' || myline[i] == '\'')
			string_starter = myline[i++];
		if (myline[i] == '\0' || myline[i] == '\n') break;
		word[wordnb++] = myline+i;	// the beginning of the word
		word[wordnb] = NULL;		// put next word to NULL
		do {
			i++;
			if (string_starter != '\0' && myline[i] == string_starter) {
				string_starter = '\0';
				break;
			}
		} while (i < len && (!isblank(myline[i]) || string_starter != '\0') && myline[i] != '\n') ;
		if (myline[i] == '\0')		// the end of the word and line (i == len)
			break;
		myline[i++] = '\0';		// the end of the word
	} while (wordnb < MAX_COMMAND_WORDS - 1) ;

	// search for the command in the list
	if (word[0] == NULL) return 1;
	i = sizeof(commande)/sizeof(command);
	while (strcasecmp (commande[--i].name, word[0])) {
		if (i == 0) {
			siril_log_message("*** Unknown command: '%s' or not implemented yet\n", word[0]);
			return 1 ;
		}
	}

	// verify argument count
	if(wordnb - 1 < commande[i].nbarg) {
		siril_log_message("   *** usage: %s\n", commande[i].usage);
		return 1;
	}

	// process the command
	commande[i].process(wordnb);

	free(myline);
	return 0;
}
/**********************************************************************/
/************** OLD FUNCTIONS: MUST BE REINTEGRATED OR NOT*************/
/**********************************************************************/

#if 0
	/* displays genname sequence from start_im to end_im with delay in 0.1 s and stride stride */
	{"animate", 3, "animate genname start_im end_im delay stride", process_animate},
	{"convert", 1, "convert filename", process_convert}, 
	/* register a sequence of images */
	{"register", 2, "register genname number [size x y]", process_register},
	/* stacking a sequence using genname.shift */
	{"composit", 2, "composit genname number shifted", process_composit},	
	{"medstack", 2, "medstack genname number outfile ", process_medstack}, 
	/* shift current image by x,y */
	{"shift", 2, "shift x y", process_shift}, 
	{"shift2", 5, "shift2 in out n x y", process_shift2}, 
	/* effectively shifts first n images sequence "in", output goes to sequence "out" */
	{"rshift2", 4, "rshift2 in out n shiftfile", process_rshift2}, 
	{"rrgb", 1, "file [file file]", process_rrgb}, /* 1 rgb file or 3 r,g,b files */
	{"grgb", 1, "file [file file]", process_grgb}, /* 1 rgb file or 3 r,g,b files */
	{"brgb", 1, "file [file file]", process_brgb}, /* 1 rgb file or 3 r,g,b files */
	{"lrgb", 2, "file file [file file]", process_lrgb},/* 1 lum file + (1 rgb file or 3 r,g,b files */
	/* i***2 commands oper the genname sequence and filename*/
	{"iadd2", 4, "iadd2 genname outname filename number", process_imoper2},
	{"isub2", 4, "isub2 genname outname filename number", process_imoper2},
	{"imul2", 4, "imul2 genname outname filename number", process_imoper2}, 
	{"idiv2", 4, "idiv2 genname outname filename number", process_imoper2},
	/* s*** commands oper scalar and curent image */
	{"sadd", 1, "sadd scalar", process_soper},
	{"smul", 1, "smul scalar", process_soper}, 
	{"sdiv", 1, "sdiv scalar", process_soper},
	{"ssub", 1, "ssub scalar", process_soper},
	/* s***2 commands oper the genname sequence and scalar */
	{"sadd2", 4, "sadd2 genname outname scalar number", process_soper2},
	{"smul2", 4, "smul2 genname outname scalar number", process_soper2}, 
	{"sdiv2", 4, "sdiv2 genname outname scalar number", process_soper2},
	{"ssub2", 4, "sdiv2 genname outname scalar number", process_soper2},
	{"composervb", 1, "composervb filename", process_composervb}, /* compose et sauve un bmp a partir de (wfit[0]), (wfit[1]), (wfit[2] */
	{"trichro", 3, "trichro rfile gfile bfile", process_trichro}, /* compose et affiche un rgb a partir de rfile, vfile, bfile */
#endif

#if 0
int process_convert(int nb){
	int type;
	char suffix[16];

	// use stat_file instead
	if (findtype(word[1], suffix, &type)){
		siril_log_message("Dont know how to convert %s.%s\n",word[1], suffix);
		return 1;
	};
	// FIXME	convert(word[1], suffix,type);
	return 0;
}

int process_trichro(int nb){
	fprintf(stderr, "process_trichro is obsolete\n");
	return 0;
}

int process_composervb(int nb){
	fprintf(stderr,"composervb: is obsolete\n");
	return 0;
}

// scalar add
int process_soper(int nb){ 
	float scalar;

	scalar=atof(word[1]);
	soper(&gfit,scalar,word[0][1]);
	redraw(com.cvport,REMAP_ALL);
	return 0;
}

int process_imoper2(int nb){
	int i;
	int number;

	number=atoi(word[4]);
	readfits(word[3], &(wfit[4]), NULL);
	for (i=1;i<=number;++i){
		readfits(buildfilename(word[1],i), &(gfit), NULL);
		process_imoper(0);
		savefits(buildfilename(word[2],i),&(gfit));
	}
	return 0;
}

int process_soper2(int nb){
	int i;
	int number;
	char in[256];

	number=atoi(word[4]);
	strncpy(in,word[1],255);
	strcpy(word[1],word[4]);
	for (i=1;i<=number;++i){
		readfits(buildfilename(in,i), &(gfit), NULL);
		process_soper(0);
		savefits(buildfilename(word[2],i),&(gfit));
	}
	return 0;
}

int process_composit(int nb){
	isempty(word[1]);
	if (nb==4){
		composit(word[1], atoi(word[2]), (char)atoi(word[3])); 
	}
	else {
		composit(word[1], atoi(word[2]), 1); 
	}
	//level_adjust(&com.g);
	redraw(com.cvport, REMAP_ALL);
	return 0;
}

int process_shift(int nb){
	shift(atoi(word[1]), atoi(word[2])); 
	redraw(com.cvport, REMAP_ALL);
	//	level_adjust(&com.g);
	return 0;
}

int process_shift2(int nb){
	int i, /*n,*/ x, y;

	isempty(word[1]);
	isempty(word[2]);
	//n=atoi(word[3]);
	x=atoi(word[4]);
	y=atoi(word[5]);

	for(i=1;i<=atoi(word[3]);++i){
		buildfilename(word[1],i);
		readfits(com.formname, &gfit, com.formname);
		shift(x,y); 
		buildfilename(word[2],i);
		savefits(com.formname,&gfit);
	}
	return 0;
}

int process_rshift2(int nb){
	isempty(word[1]);
	isempty(word[2]);
	if (nb==5){
		rshift2(word[1], word[2], atoi(word[3]), word[4]); 
	}
	level_adjust(&gfit);
	return 0;
}

int process_unsharp2(int nb){
	int i;

	for (i=0;i<=atoi(word[5]);++i){
		if(!readfits(buildfilename(word[3], i), &(gfit), NULL)){
			process_unsharp(nb);
			savefits(buildfilename(word[4], i), &(gfit));
		}
	}
	return 0;
}

int process_crop2(int nb){
	int i, err;
	char ofilename[256];

	if (!com.drawn){
		if (nb==8){
			com.rectX = atoi(word[4]);
			com.rectY = atoi(word[5]);
			com.rectW = atoi(word[6]);
			com.rectH = atoi(word[7]);
		}
		else {
			siril_log_message("Crop2: select a region or provide x,y,width,height\n");
			return 1;
		}
	}

	for (i=1;i<=atoi(word[3]);++i){	
		buildfilename(word[1],i);
		snprintf(ofilename,255,"%s%d",word[2],i);
		err=readfits(com.formname,&(gfit), NULL);
		if (err){
			siril_log_message("Crop2: aborted\n");
			return 1;
		}
		crop(&(gfit));
		savefits(ofilename,&(gfit));
	}
	redraw(com.cvport, REMAP_ALL);
	return 0;
}

int process_medstack(int nb){
	char filename[256];
	stat_file(word[3], &(com.imagetype), filename);
	medstack(word[1],atoi(word[2]),filename);
	strcpy(word[1],filename);
	process_load(1);
	return 0;
}

int process_register(int nb){
	struct registration_args reg_args;

	if (com.drawn){
		reg_args.sel1size = max(com.rectW, com.rectH);
		reg_args.sel1X = com.rectX + com.rectW/2;
		reg_args.sel1Y = com.rectY + com.rectH/2;
	}
	else{
		if (nb<6){
			siril_log_message("Select a rectangle or provide centerx centery size\n");
			return 0;
		}
		reg_args.sel1X = atoi(word[3]);
		reg_args.sel1Y = atoi(word[4]);
		reg_args.sel1size = atoi(word[5]);
	}
	reg_args.all_images = TRUE;
	//reg_args.number = atoi(word[2]);
	isempty(word[1]);
	readseqfile(com.seq.name);
	register_shift(word[1], &reg_args, com.reglayer); 
	writeseqfile(com.seq.name);
	return 0;
}

int process_animate(int nb){
	int i, start, end, stride, err;
	double delay;

	start=atoi(word[2]);
	end=atoi(word[3]);
	delay=atof(word[4]);
	stride=max(1,atoi(word[5]));
	if(delay <0.02)
		delay=0.02;
	if(delay >1)
		delay=(double)1;
	for (i=start;i<=end && !(com.busy & ST_CANCEL);i+=stride){
		timing(0,"");
		buildfilename(word[1],i);
		err=readfits(com.formname, &(gfit), com.formname);
		siril_log_message("Displaying image %s %f\n", com.formname,delay);
		progress(0);
		if (err){
			siril_log_message("animate: error\n");
		}
		else{
			redraw(com.cvport,REMAP_ALL);
		}

		while (timing(1,"")<delay){
			progress(0);
		}
	}
	return 0;
}

int xrgb(int nb, int layer){
	siril_log_message("Entering xrgb command, nb %d reading %s:%s:%s\n", nb, word[1], word[2], word[3]);
	switch(nb){
		case 2: // only one RGB filename is given
			// according to the value x of layer
			// perform xrgb composite
			readfits(word[1],wfit+4, NULL);
			copyfits(wfit+4,&gfit,CP_ALLOC|CP_EXPAND|CP_FORMAT,0); 	/* lrgb result goes to gfit */
			copyfits(wfit+4,wfit,CP_ALLOC|CP_EXTRACT,layer); 	/* l */
			copyfits(wfit+4,wfit+1,CP_ALLOC|CP_EXTRACT,RLAYER);	/* r */
			copyfits(wfit+4,wfit+2,CP_ALLOC|CP_EXTRACT,GLAYER);	/* g */
			copyfits(wfit+4,wfit+3,CP_ALLOC|CP_EXTRACT,BLAYER);	/* b */
			break;

		case 4:
			readfits(word[1],wfit+1,NULL);	/* r */
			copyfits(wfit+1,&gfit,CP_ALLOC|CP_EXPAND|CP_FORMAT,0); 	/* lrgb result goes to gfit */
			readfits(word[2],wfit+2,NULL);	/* g */
			readfits(word[3],wfit+3,NULL);	/* b */
			copyfits(wfit+layer+1,wfit,CP_COPYA|CP_ALLOC|CP_FORMAT,layer); 	/* l */
			break;

		default:
			// error
			break;
	}
	fprintf(stderr,"invoking lrgb command from xrgb\n");
	lrgb(wfit, wfit+1, wfit+2, wfit+3, &gfit);
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_rrgb(int nb){
	fprintf(stderr,"entering rrgb command\n");
	xrgb(nb,RLAYER);
	return 0;
}

int process_grgb(int nb){
	fprintf(stderr,"entering grgb command\n");
	xrgb(nb,GLAYER);
	return 0;
}

int process_brgb(int nb){
	fprintf(stderr,"entering brgb command\n");
	xrgb(nb,BLAYER);
	return 0;
}

int process_lrgb(int nb){
	//
	// lrgb combines an rgb fits and a luminance fits
	// through a short irruption in HSI space.
	// takes 2 or 3 arguments. 
	//			lrgb with 2 filenames l,rgb as arg1 and arg2
	//			lrgb with 4 filenames l,r,g,b as arg[1234]
	// 	
	fprintf(stderr,"entering lrgb command\n");
	switch(nb){
		case 2:
			readfits(word[2],wfit,NULL);
			readfits(word[3],wfit+4,NULL);
			copyfits(wfit+4,wfit+1,CP_EXTRACT,RLAYER); 	/* r */
			copyfits(wfit+4,wfit+2,CP_EXTRACT,GLAYER);	/* g */
			copyfits(wfit+4,wfit+3,CP_EXTRACT,BLAYER);	/* b */
			break;

		case 3:
			readfits(word[2],wfit,NULL);
			readfits(word[3],wfit+4,NULL);
			readfits(word[4],wfit+4,NULL);
			readfits(word[5],wfit+4,NULL);
			break;	

		default:
			// error
			break;
	}
	fprintf(stderr,"invoking lrgb function\n");
	copyfits(wfit,&gfit,CP_ALLOC|CP_EXPAND|CP_FORMAT,0); 	/* lrgb result goes to gfit */

	lrgb(wfit, wfit+1, wfit+2, wfit+3, &gfit);	
	redraw(com.cvport, REMAP_ALL);
	return 0;
}
#endif
