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
#include <config.h>
#endif

#include <gtk/gtk.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/time.h>
#include <dirent.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <libgen.h>

#include "siril.h"
#include "proto.h"
#include "callbacks.h"
#include "ser.h"
#ifdef HAVE_FFMS2
#include "films.h"
#endif
#include "single_image.h"
#include "histogram.h"
#include "PSF.h"
#include "PSF_list.h"	// clear_stars_list
#include "quality.h"
#include "registration.h"	// for update_reg_interface
#include "stacking.h"	// for update_stack_interface

/* when opening a file outside the main sequence loading system and that file
 * is a sequence (SER/AVI), this function is called to load this sequence. */
int read_single_sequence(char *realname, int imagetype) {
	int retval=3;		// needs to return 3 if ok !!!
	char *name = strdup(realname);
	char *dirname = extract_path(realname);
	if (!changedir(dirname))
		writeinitfile();
	free(dirname);

	if (check_only_one_film_seq(realname)) retval = 1;
	else {
#ifdef HAVE_FFMS2
	const char *ext;
#endif
		switch (imagetype) {
			case TYPESER:
				name[strlen(name)-1] = 'q';
				break;
#ifdef HAVE_FFMS2
			case TYPEAVI:
				ext = get_filename_ext(realname);
				int len = strlen(ext);
				strncpy(name+strlen(name)-len, "seq", len);
				break;
#endif
			default:
				retval = 1;
		}
		if (!set_seq(basename(name))) {
			/* if it loads, make it selected and only element in the list of sequences */
			control_window_switch_to_tab(IMAGE_SEQ);
			GtkComboBoxText *combo_box_text = GTK_COMBO_BOX_TEXT(lookup_widget("sequence_list_combobox"));
			gtk_combo_box_text_remove_all(combo_box_text);
			gtk_combo_box_text_append(combo_box_text, 0, basename(realname));
			g_signal_handlers_block_by_func(GTK_COMBO_BOX(combo_box_text), on_seqproc_entry_changed, NULL);
			gtk_combo_box_set_active(GTK_COMBO_BOX(combo_box_text), 0);
			g_signal_handlers_unblock_by_func(GTK_COMBO_BOX(combo_box_text), on_seqproc_entry_changed, NULL);
		}
		else retval = 1;
	}
	free(name);
	return retval;
}

/* Find sequences in CWD and create .seq files.
 * In the current working directory, looks for sequences of fits files or files
 * already representing sequences like SER and AVI formats and builds the
 * corresponding sequence files.
 * Called when changing wd with name == NULL or when an explicit root name is
 * given in the GUI or when searching for sequences.
 */
int check_seq(int force) {
	char *basename;
	int curidx, fixed;
	DIR *dir;
	struct dirent *file;
	sequence *sequences[40];
	int i, nb_seq = 0;

	if (!com.wd) {
		siril_log_message("Current working directory is not set, aborting.\n");
		return 1;
	}
	if ((dir = opendir(com.wd)) == NULL) {
		fprintf(stderr, "working directory cannot be opened.\n");
		free(com.wd);
		com.wd = NULL;
		return 1;
	}

	while ((file = readdir(dir)) != NULL) {
		sequence *new_seq;
		int fnlen = strlen(file->d_name);
		if (!strncasecmp(file->d_name + fnlen - 4, ".ser", 4)) {
			struct ser_struct *ser_file = malloc(sizeof(struct ser_struct));
			ser_init_struct(ser_file);
			if (ser_open_file(file->d_name, ser_file))
				continue;
			new_seq = calloc(1, sizeof(sequence));
			initialize_sequence(new_seq, TRUE);
			new_seq->seqname = strndup(file->d_name, fnlen-4);
			new_seq->beg = 0;
			new_seq->end = ser_file->frame_count-1;
			new_seq->number = ser_file->frame_count;
			new_seq->type = SEQ_SER;
			new_seq->ser_file = ser_file;
			sequences[nb_seq] = new_seq;
			nb_seq++;
			fprintf(stdout, "Found a SER sequence (number %d)\n", nb_seq);
			set_progress_bar_data(NULL, PROGRESS_PULSATE);
		}
#ifdef HAVE_FFMS2
		else if (!check_for_film_extensions(get_filename_ext(file->d_name))) {
			struct film_struct *film_file = malloc(sizeof(struct film_struct));
			if (film_open_file(file->d_name, film_file)) {
				free(film_file);
				continue;
			}
			new_seq = calloc(1, sizeof(sequence));
			initialize_sequence(new_seq, TRUE);
			int len = strlen(get_filename_ext(file->d_name));
			new_seq->seqname = strndup(file->d_name, fnlen-(len+1));
			new_seq->beg = 0;
			new_seq->end = film_file->frame_count-1;
			new_seq->number = film_file->frame_count;
			new_seq->type = SEQ_AVI;
			new_seq->film_file = film_file;
			sequences[nb_seq] = new_seq;
			nb_seq++;
			fprintf(stdout, "Found a AVI sequence (number %d)\n", nb_seq);
			set_progress_bar_data(NULL, PROGRESS_PULSATE);
		}
#endif

		else if (!strncasecmp(file->d_name + fnlen - com.len_ext, com.ext,
				com.len_ext)) {
			if (!get_index_and_basename(file->d_name, &basename, &curidx, &fixed)) {
				int current_seq = -1;
				/* search in known sequences if we already have it */
				for (i=0; i<nb_seq; i++) {
					if (!strcmp(sequences[i]->seqname, basename)) {
						current_seq = i;
					}
				}
				/* not found */
				if (current_seq == -1) {
					if (nb_seq == 40) {
						fprintf(stderr, "too many sequences\n");
						continue;
					}
					new_seq = calloc(1, sizeof(sequence));
					initialize_sequence(new_seq, TRUE);
					new_seq->seqname = basename;
					new_seq->beg = INT_MAX;
					new_seq->end = 0;
					new_seq->fixed = fixed;
					sequences[nb_seq] = new_seq;
					current_seq = nb_seq;
					nb_seq++;
					fprintf(stdout, "Found a sequence (number %d) with base name"
							" \"%s\", looking for first and last indexes.\n",
							nb_seq, basename);
					set_progress_bar_data(NULL, PROGRESS_PULSATE);
				}
				if (curidx < sequences[current_seq]->beg)
					sequences[current_seq]->beg = curidx;
				if (curidx > sequences[current_seq]->end)
					sequences[current_seq]->end = curidx;
			}
		}
	}
	closedir(dir);
	if (nb_seq > 0) {
		int retval = 1;
		for (i=0; i<nb_seq; i++) {
			if (sequences[i]->beg != sequences[i]->end) {
				char msg[200];
				sprintf(msg, "sequence %d, found: %d to %d",
						i+1, sequences[i]->beg, sequences[i]->end);
				set_progress_bar_data(msg, PROGRESS_NONE);
				if (!buildseqfile(sequences[i], force) && retval)
					retval = 0;	// at least one succeeded to be created
			}
			free_sequence(sequences[i], TRUE);
		}
		return retval;
	}
	return 1;	// no sequence found
}

/* Check for on film sequence of the name passed in arguement
 * Returns 0 if OK */
int check_only_one_film_seq(char* name) {
	int retval = 1;
	DIR *dir;
	sequence *new_seq = NULL;

	if (!com.wd) {
		siril_log_message("Current working directory is not set, aborting.\n");
		return 1;
	}
	if ((dir = opendir(com.wd)) == NULL) {
		fprintf(stderr, "working directory cannot be opened.\n");
		free(com.wd);
		com.wd = NULL;
		return 1;
	}
	
	int fnlen = strlen(name);
	if (!strncasecmp(name + fnlen - 4, ".ser", 4)) {
		struct ser_struct *ser_file = malloc(sizeof(struct ser_struct));
		ser_init_struct(ser_file);
		if (ser_open_file(name, ser_file)) return 1;
			
		new_seq = calloc(1, sizeof(sequence));
		initialize_sequence(new_seq, TRUE);
		new_seq->seqname = strndup(name, fnlen-4);
		new_seq->beg = 0;
		new_seq->end = ser_file->frame_count-1;
		new_seq->number = ser_file->frame_count;
		new_seq->type = SEQ_SER;
		new_seq->ser_file = ser_file;
	}
#ifdef HAVE_FFMS2
	else if (!check_for_film_extensions(get_filename_ext(name))) {
		struct film_struct *film_file = malloc(sizeof(struct film_struct));
		if (film_open_file(name, film_file)) {
			free(film_file);
			return 1;
		}
		new_seq = calloc(1, sizeof(sequence));
		initialize_sequence(new_seq, TRUE);
		int len = strlen(get_filename_ext(name));
		new_seq->seqname = strndup(name, fnlen-len-1);
		new_seq->beg = 0;
		new_seq->end = film_file->frame_count-1;
		new_seq->number = film_file->frame_count;
		new_seq->type = SEQ_AVI;
		new_seq->film_file = film_file;
		fprintf(stdout, "Found a AVI sequence\n");
	}
#endif
	closedir(dir);
	if (!new_seq) return 0;
	if (new_seq->beg != new_seq->end) {
			if (!buildseqfile(new_seq, 0) && retval)
				retval = 0;
		}
		free_sequence(new_seq, TRUE);
	return retval;
}

/* load a sequence and initialized everything that relates */
int set_seq(const char *name){
	sequence *seq;
	int image_to_load;
	
	if ((seq = readseqfile(name)) == NULL) {
		fprintf(stderr, "could not load sequence %s\n", name);
		return 1;
	}
	free_image_data();
	if (seq->reference_image != -1)
		image_to_load = seq->reference_image;
	else image_to_load = 0;

	if (seq_read_frame(seq, image_to_load, &gfit)) {
		fprintf(stderr, "could not load first image from sequence\n");
		free(seq);
		return 1;
	}

	/* initialize sequence-related runtime data */
	seq->rx = gfit.rx; seq->ry = gfit.ry;
	seq->current = image_to_load;

	if (seq->nb_layers == -1 || seq->nb_layers != gfit.naxes[2]) {	// not init yet, first loading of the sequence
		seq->nb_layers = gfit.naxes[2];
		seq->regparam = calloc(seq->nb_layers, sizeof(regdata *));
		seq->layers = calloc(seq->nb_layers, sizeof(layer_info));
		writeseqfile(seq);
	}
	/* Sequence is stored in com.seq for now */
	free_sequence(&com.seq, FALSE);
	memcpy(&com.seq, seq, sizeof(sequence));

	if (seq->nb_layers > 1)
		show_rgb_window();
	else hide_rgb_window();
	init_layers_hi_and_lo_values(MIPSLOHI); // set some hi and lo values in seq->layers,
	set_cutoff_sliders_max_values();// update min and max values for contrast sliders
	set_cutoff_sliders_values();	// update values for contrast sliders for this image
	seqsetnum(image_to_load);	// set limits for spin button and display loaded filenum
	set_layers_for_assign();	// set default layers assign and populate combo box
	set_layers_for_registration();	// set layers in the combo box for registration
	fill_sequence_list(seq, 0);	// display list of files in the sequence
	set_output_filename_to_sequence_name();
	sliders_mode_set_state(com.sliders);
	initialize_display_mode();

	/* initialize image-related runtime data */
	set_display_mode();		// display the display mode in the combo box
	display_filename();		// display filename in gray window
	adjust_exclude(image_to_load, FALSE);	// check or uncheck excluded checkbox
	adjust_refimage(image_to_load);	// check or uncheck reference image checkbox
	set_prepro_button_sensitiveness(); // enable or not the preprobutton
	update_reg_interface(TRUE);	// change the registration prereq message
	update_stack_interface();	// get stacking info and enable the Go button
	adjust_reginfo();		// change registration displayed/editable values
	update_gfit_histogram_if_needed();
	adjust_sellabel();

	/* redraw and display image */
	show_main_gray_window();
	close_tab();	//close Green and Blue Tab if a 1-layer sequence is loaded
	adjust_vport_size_to_image();	// resize viewports to the displayed image size
	redraw(com.cvport, REMAP_ALL);

	update_used_memory();
	return 0;
}

/* Load image number index from the sequence and display it.
 * if load_it is true, dest is assumed to be gfit
 * TODO: cut that method in two, with an internal func taking a filename and a fits
 */
int seq_load_image(sequence *seq, int index, fits *dest, gboolean load_it) {
	seq->current = index;
	clear_stars_list();
	clear_histograms();
	gfit.maxi = 0;
	// what else needs to be cleaned?
	if (load_it) {
		set_cursor_waiting(TRUE);
		if (seq_read_frame(seq, index, dest)) {
			set_cursor_waiting(FALSE);
			return 1;
		}
		set_fwhm_star_as_star_list(seq);// display the fwhm star if possible
		if (com.sliders != USER) {
			init_layers_hi_and_lo_values(com.sliders);
			sliders_mode_set_state(com.sliders);
			set_cutoff_sliders_max_values();// update min and max values for contrast sliders
			set_cutoff_sliders_values();	// update values for contrast sliders for this image
			set_display_mode();		// display the display mode in the combo box
		}
		redraw(com.cvport, REMAP_ALL);	// redraw and display image
		redraw_previews();		// redraw registration preview areas
		display_filename();		// display filename in gray window
		adjust_reginfo();		// change registration displayed/editable values
		calculate_fwhm(com.vport[com.cvport]);
		update_gfit_histogram_if_needed();
		set_cursor_waiting(FALSE);
	}
	/* change the displayed value in the spin button to have the real file number
	 * instead of the index of the adjustment */
	display_image_number(index);
	sequence_list_change_current();
	adjust_exclude(index, FALSE);	// check or uncheck excluded checkbox
	adjust_refimage(index);	// check or uncheck reference image checkbox
	update_used_memory();
	return 0;
}

/*****************************************************************************
 *              SEQUENCE FUNCTIONS FOR NON-OPENED SEQUENCES                  *
 * **************************************************************************/

/* Get the filename of an image in a sequence.
 * Return value is the same as the name_buf argument, which must be
 * pre-allocated to at least 256 characters. If sequence has no file names, a
 * description like image "42 from awesome_mars.ser" is made. */
char *seq_get_image_filename(sequence *seq, int index, char *name_buf) {
	switch (seq->type) {
		case SEQ_REGULAR:
			return fit_sequence_get_image_filename(seq, index, name_buf, TRUE);
		case SEQ_SER:
			if (!name_buf || index < 0 || index > seq->end) {
				return NULL;
			}
			snprintf(name_buf, 255, "%d from %s.ser", index, seq->seqname);
			name_buf[255] = '\0';
			return name_buf;
#ifdef HAVE_FFMS2
		case SEQ_AVI:
			if (!name_buf || index < 0 || index > seq->end) {
				return NULL;
			}
			//snprintf(name_buf, 255, "%d from %s.avi", index, seq->seqname);
			snprintf(name_buf, 255, "%s_%d", seq->seqname, index);
			name_buf[255] = '\0';
			return name_buf;
#endif
		case SEQ_INTERNAL:
			snprintf(name_buf, 255, "%s_%d", seq->seqname, index);
			name_buf[255] = '\0';
			return name_buf;
	}
	return NULL;
}

/* Read an entire image from a sequence, inside a pre-allocated fits.
 * Opens the file, reads data, closes the file.
 */
int seq_read_frame(sequence *seq, int index, fits *dest) {
	char filename[256];
	assert(index < seq->number);
	switch (seq->type) {
		case SEQ_REGULAR:
			fit_sequence_get_image_filename(seq, index, filename, TRUE);
			if (readfits(filename, dest, NULL)) {
				siril_log_message("could not load image %d from sequence %s\n",
						index, seq->seqname); 
				return 1;
			}
			break;
		case SEQ_SER:
			assert(seq->ser_file);
			if (ser_read_frame(seq->ser_file, index, dest)) {
				siril_log_message("could not load frame %d from SER sequence %s\n",
						index, seq->seqname); 
				return 1;
			}
			break;
#ifdef HAVE_FFMS2
		case SEQ_AVI:
			assert(seq->film_file);
			if (film_read_frame(seq->film_file, index, dest)) {
				siril_log_message("could not load frame %d from AVI sequence %s\n",
						index, seq->seqname); 
				return 1;
			}
			// should dest->maxi be set to 255 here?
			break;
#endif
		case SEQ_INTERNAL:
			assert(seq->internal_fits);
			copyfits(seq->internal_fits[index], dest, CP_FORMAT, -1);
			dest->data = seq->internal_fits[index]->data;
			dest->pdata[0] = seq->internal_fits[index]->pdata[0];
			dest->pdata[1] = seq->internal_fits[index]->pdata[1];
			dest->pdata[2] = seq->internal_fits[index]->pdata[2];
			break;
	}
	image_find_minmax(dest, 0);
	return 0;
}

/* same as seq_read_frame above, but creates an image the size of the selection
 * rectangle only. layer is set to the layer number in the read partial frame.
 * The partial image result is only one-channel deep, so it cannot be used to
 * have a partial RGB image. */
int seq_read_frame_part(sequence *seq, int layer, int index, fits *dest, const rectangle *area) {
	char filename[256];
	fits tmp_fit;
	memset(&tmp_fit, 0, sizeof(fits));
	switch (seq->type) {
		case SEQ_REGULAR:
			fit_sequence_get_image_filename(seq, index, filename, TRUE);
			if (readfits_partial(filename, layer, dest, area)) {
				siril_log_message("Could not load partial image %d from sequence %s\n",
						index, seq->seqname); 
				return 1;
			}
			break;
		case SEQ_SER:
			assert(seq->ser_file);
			/* TODO: build a FITS from ser_read_opened_partial() */
			if (ser_read_frame(seq->ser_file, index, &tmp_fit)) {
				siril_log_message("could not load frame %d from SER sequence %s\n",
						index, seq->seqname); 
				return 1;
			}
			extract_region_from_fits(&tmp_fit, layer, dest, area);
			clearfits(&tmp_fit);
			break;
#ifdef HAVE_FFMS2
		case SEQ_AVI:
			assert(seq->film_file);
			if (film_read_frame(seq->film_file, index, &tmp_fit)) {
				siril_log_message("could not load frame %d from AVI sequence %s\n",
						index, seq->seqname); 
				return 1;
			}
			extract_region_from_fits(&tmp_fit, layer, dest, area);
			clearfits(&tmp_fit);
			break;
#endif
		case SEQ_INTERNAL:
			assert(seq->internal_fits);
			extract_region_from_fits(seq->internal_fits[index], 0, dest, area);
			break;
	}
	return 0;
}

/*****************************************************************************
 *                 SEQUENCE FUNCTIONS FOR OPENED SEQUENCES                   *
 * **************************************************************************/

/* locks cannot be probed to see if they are init or not, so we have to keep
 * all of them in the same state, which is initialized if the array is non-nul. */
static int _allocate_sequence_locks(sequence *seq) {
#ifdef _OPENMP
	if (!seq->fd_lock) {
		int i;
		seq->fd_lock = malloc(seq->number * sizeof(omp_lock_t));
		if (!seq->fd_lock) {
			fprintf(stderr, "Allocation error when opening images, aborting\n");
			return 1;
		}

		for (i=0; i<seq->number; i++)
			omp_init_lock(&seq->fd_lock[i]);
	}
#endif
	return 0;
}

/* open image for future intensive operations (read only) */
int seq_open_image(sequence *seq, int index) {
	int status = 0;
	char filename[256];
	switch (seq->type) {
		case SEQ_REGULAR:
			if (!seq->fptr) {
				seq->fptr = calloc(seq->number, sizeof(fitsfile *));
				if (!seq->fptr) {
				       fprintf(stderr, "Allocation error when opening images, aborting\n");
			       	       return 1;
				}
			}
			if (_allocate_sequence_locks(seq))
				return 1;

			fit_sequence_get_image_filename(seq, index, filename, TRUE);
			fits_open_file(&seq->fptr[index], filename, READONLY, &status);
			if (status) {
				fits_report_error(stderr, status);
				return status;
			}
			/* should we check image parameters here? such as bitpix or naxis */
			break;
		case SEQ_SER:
			assert(seq->ser_file->fd > 0);
			break;
#ifdef HAVE_FFMS2
		case SEQ_AVI:
			siril_log_message("This operation is not supported on AVI sequences (seq_open_image)\n");
			return 1;
#endif
		case SEQ_INTERNAL:
			siril_log_message("This operation is not supported on internal sequences (seq_open_image)\n");
			return 1;
	}
	return 0;
}

/* close opened images, only useful for regular FITS sequences */
void seq_close_image(sequence *seq, int index) {
	int status = 0;
	switch (seq->type) {
		case SEQ_REGULAR:
			if (seq->fptr && seq->fptr[index]) {
				fits_close_file(seq->fptr[index], &status);
				seq->fptr[index] = NULL;
			}
			break;
		default:
			break;
	}
}

/* read a region in a layer of an opened file from a sequence.
 * The buffer must have been allocated to the size of the area. */
int seq_opened_read_region(sequence *seq, int layer, int index, WORD *buffer, const rectangle *area) {
	switch (seq->type) {
		case SEQ_REGULAR:
			return read_opened_fits_partial(seq, layer, index, buffer, area);
		case SEQ_SER:
			return ser_read_opened_partial(seq->ser_file, layer, index, buffer, area);
		default:
			break;
	}
	return 0;
}


/*****************************************************************************
 *                         SEQUENCE DATA MANAGEMENT                          *
 * **************************************************************************/

/* if FWHM was calculated on the sequence, a minimisation exists for all
 * images, and when switching to a new image, it should be set as the only item
 * in the star list, in order to be displayed.
 * A special care is required in PSF_list.c:clear_stars_list(), to not free this data. */
void set_fwhm_star_as_star_list_with_layer(sequence *seq, int layer) {
	assert(seq->regparam);
	/* we chose here the first layer that has been allocated, which doesn't
	 * mean it contains data for all images. Handle with care. */
	if (seq->regparam && layer >= 0 && layer < seq->nb_layers && seq->regparam[layer] &&
			seq->regparam[layer][seq->current].fwhm_data && !com.stars) {
		com.stars = malloc(2 * sizeof(fitted_PSF *));
		com.stars[0] = seq->regparam[layer][seq->current].fwhm_data;
		com.stars[1] = NULL;
		com.star_is_seqdata = TRUE;
	}
}

// cannot be called in the worker thread
void set_fwhm_star_as_star_list(sequence *seq) {
	int layer = get_registration_layer(seq);
	set_fwhm_star_as_star_list_with_layer(seq, layer);
}

/* Rebuilds the file name of an image in a sequence.
 * The file name is stored in name_buffer, which must be allocated 256 bytes
 * The index is the index in the sequence, not the number appearing in the file name
 * Return value: NULL on error, name_buffer on success.
 */
char *fit_sequence_get_image_filename(sequence *seq, int index, char *name_buffer, gboolean add_fits_ext) {
	char format[20];
	if (index < 0 || index > seq->number || name_buffer == NULL)
		return NULL;
	if (seq->fixed <= 1){
		sprintf(format, "%%s%%d");
	} else {
		sprintf(format, "%%s%%.%dd", seq->fixed);
	}
	if (add_fits_ext)
		strcat(format, com.ext);
	snprintf(name_buffer, 255, format,
			seq->seqname, seq->imgparam[index].filenum);
	name_buffer[255] = '\0';
	return name_buffer;
}

/* Returns a filename for an image that could be in a sequence, but the sequence structure
 * has not been fully initialized yet. Only beg, end, fixed and seqname are used.
 */
char *get_possible_image_filename(sequence *seq, int image_number, char *name_buffer) {
	char format[20];
	if (image_number < seq->beg || image_number > seq->end || name_buffer == NULL)
		return NULL;
	if (seq->fixed <= 1){
		sprintf(format, "%%s%%d%s", com.ext);
	} else {
		sprintf(format, "%%s%%.%dd%s", seq->fixed, com.ext);
	}
	sprintf(name_buffer, format, seq->seqname, image_number);
	return name_buffer;
}

/* splits a filename in a base name and an index number, if the file name ends with .fit
 * it also computes the fixed length if there are zeros in the index */
int	get_index_and_basename(const char *filename, char **basename, int *index, int *fixed){
	char *buffer;
	int i, fnlen, first_zero, digit_idx;

	*index = -1;		// error values
	*fixed = 0;
	first_zero = -1;
	*basename = NULL;
	fnlen = strlen(filename);
	if (fnlen < strlen(com.ext)+2) return -1;
	if (!ends_with(filename, com.ext)) return -1;
	i = fnlen-strlen(com.ext)-1;
	if (!isdigit(filename[i])) return -1;
	digit_idx = i;

	buffer = strdup(filename);
	buffer[fnlen-strlen(com.ext)] = '\0';		// for atoi()
	do {
		if (buffer[i] == '0' && first_zero < 0)
			first_zero = i;
		if (buffer[i] != '0' && first_zero > 0)
			first_zero = -1;
		i--;
	} while (i >= 0 && isdigit(buffer[i]));
	i++;
	if (i == 0) {
		free(buffer);
		return -1;	// no base name, only number
	}
	if (first_zero >= 0)
		*fixed = digit_idx - i + 1;
	//else *fixed = 0;
	*index = atoi(buffer+i);
	if (*basename == NULL) {	// don't copy it if we already have it
		*basename = malloc(i * sizeof(char) + 1);
		strncpy(*basename, buffer, i);
		(*basename)[i] = '\0';
	}
	//fprintf(stdout, "from filename %s, base name is %s, index is %d\n", filename, *basename, *index);
	free(buffer);
	return 0;
}

/* sets default values for the sequence */
void initialize_sequence(sequence *seq, gboolean is_zeroed) {
	int i;
	if (!is_zeroed) {
		memset(seq, 0, sizeof(sequence));
	}
	seq->nb_layers = -1;		// uninit value
	seq->reference_image = -1;	// uninit value
	seq->type = SEQ_REGULAR;
	for (i=0; i<PREVIEW_NB; i++) {
		seq->previewX[i] = -1;
		seq->previewY[i] = -1;
	}
}

/* call this to close a sequence. Second arg must be FALSE for com.seq
 * WARNING: the data is not reset to NULL, if seq is to be reused,
 * initialize_sequence() must be called on it right after free_sequence()
 * (= do it for com.seq) */
void free_sequence(sequence *seq, gboolean free_seq_too) {
	static GtkComboBoxText *cbbt_layers = NULL;
	int i;
		
	if (cbbt_layers == NULL)
		cbbt_layers = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(
					builder, "comboboxreglayer"));
	gtk_combo_box_text_remove_all(cbbt_layers);
	
	if (seq == NULL) return;
	if (seq->nb_layers > 0 && seq->regparam) {
		for (i=0; i<seq->nb_layers; i++) {
			if (seq->regparam[i]) {
				int j;
				for (j=0; j < seq->number; j++) {
					if (seq->regparam[i][j].fwhm_data)
						free(seq->regparam[i][j].fwhm_data);
				}
				free(seq->regparam[i]);
			}
		}
		free(seq->regparam);
	}

	for (i=0; i<seq->number; i++) {
		if (seq->fptr && seq->fptr[i]) {
			int status = 0;
			fits_close_file(seq->fptr[i], &status);
		}
		if (seq->imgparam && seq->imgparam[i].stats) {
			free(seq->imgparam[i].stats);
		}
	}
	if (seq->seqname)	free(seq->seqname);
	if (seq->layers)	free(seq->layers);
	if (seq->imgparam)	free(seq->imgparam);
	if (seq->fptr)		free(seq->fptr);

#ifdef _OPENMP
	if (seq->fd_lock) {
		for (i=0; i<seq->number; i++) {
			omp_destroy_lock(&seq->fd_lock[i]);
		}
		free(seq->fd_lock);
	}
#endif

	if (seq->ser_file) {
		ser_close_file(seq->ser_file);	// frees the data too
		free(seq->ser_file);
	}
#ifdef HAVE_FFMS2
	if (seq->film_file) {
		film_close_file(seq->film_file);	// frees the data too
		free(seq->film_file);
	}
#endif
	if (seq->internal_fits) {
		/* the fits in internal_fits should still be referenced somewhere */
		free(seq->internal_fits);
	}
	if (free_seq_too)	free(seq);
}

void sequence_free_preprocessing_data(sequence *seq) {
	// free opened files
	if (seq->ppprefix) {
		free(seq->ppprefix);
		seq->ppprefix = NULL;
	}
	if (seq->offset) {
		clearfits(seq->offset);
		free(seq->offset);
		seq->offset = NULL;
	}
	if (seq->dark) {
		clearfits(seq->dark);
		free(seq->dark);
		seq->dark = NULL;
	}
	if (seq->flat) {
		clearfits(seq->flat);
		free(seq->flat);
		seq->flat = NULL;
	}
}

gboolean sequence_is_loaded() {
	return (com.seq.seqname != NULL && com.seq.imgparam != NULL);
}

/*****************************************************************************
 *                             SEQUENCE PROCESSING                           *
 * **************************************************************************/

/* Start a processing on all images of the sequence seq, on layer layer if it applies.
 * The see coment in siril.h for help on process format.
 */
int sequence_processing(sequence *seq, sequence_proc process, int layer) {
	int i;
	float cur_nb, nb_frames;
	fits fit;
	rectangle area;

	if (!com.selection.w || !com.selection.h) {
		siril_log_message("No selection was made for a selection-based sequence processing\n");
		return 1;
	}
	memcpy(&area, &com.selection, sizeof(rectangle));
	memset(&fit, 0, sizeof(fits));
	check_or_allocate_regparam(seq, layer);

	/*if (selected_images_only)
		nb_frames = (float)seq->selnum;
	else*/
	nb_frames = (float)seq->number;

	for (i=0, cur_nb=0.f; i<seq->number; ++i) {
		if (!get_thread_run()) break;
		//if (selected_images_only && !seq->imgparam[i].incl)
		//	continue;
		/* opening the image */
		if (seq_read_frame_part(seq, layer, i, &fit, &area))
			return 1;
		/* processing the image */
		if (process(seq, layer, i, &fit) < 0)
			return 1;
		cur_nb += 1.f;
		set_progress_bar_data(NULL, cur_nb/nb_frames);
	}
	return 0;
}

/* Computes FWHM for a sequence image and store data in the sequence imgdata.
 * seq_layer is the corresponding layer in the raw image from the sequence.
 */
int seqprocess_fwhm(sequence *seq, int seq_layer, int frame_no, fits *fit) {
	rectangle area;
	area.x = area.y = 0;
	area.w = fit->rx; area.h = fit->ry;
	assert(seq_layer < seq->nb_layers);
	fitted_PSF *result = get_Minimisation(fit, 0, &area);
	if (result) {
		seq->regparam[seq_layer][frame_no].fwhm_data = result;
		seq->regparam[seq_layer][frame_no].fwhm = result->fwhmx;
		assert(result->fwhmx >= result->fwhmy);
		/* just a try, to verify what's written elsewhere. If it fails,
		 * we need to add a test and store FWHMY in fwhm instead. */
		return 0;
	} else {
		seq->regparam[seq_layer][frame_no].fwhm_data = NULL;
		seq->regparam[seq_layer][frame_no].fwhm = 0.0f;
		return 1;
	}
}

void do_fwhm_sequence_processing(sequence *seq, int layer) {
	int i, retval;
	siril_log_message("Starting sequence processing of PSF\n");
	set_progress_bar_data("Computing PSF on selected star", PROGRESS_NONE);
	retval = sequence_processing(seq, &seqprocess_fwhm, layer);
	siril_log_message("Finished sequence processing of PSF\n");
	if (retval) {
		set_progress_bar_data("Failed to compute PSF for the sequence. Ready.", PROGRESS_NONE);
		set_cursor_waiting(FALSE);
		return;
	}
	// update the list
	if (seq->type != SEQ_INTERNAL)
		fill_sequence_list(seq, layer);
	// set the coordinates of the detected star, not known by the processing
	for (i=0; i<seq->number; i++) {
		fitted_PSF *star = seq->regparam[layer][i].fwhm_data;
		if (star) {
			// same code as in add_star(), set position of star using selection coords
			star->xpos = star->x0 + com.selection.x;
			star->ypos = com.selection.y + com.selection.h - star->y0;
			fprintf(stdout, "star image %d: %g, %g\n", i, star->xpos, star->ypos);
		}
	}
	set_fwhm_star_as_star_list_with_layer(seq, layer);
	set_progress_bar_data("Finished computing PSF for the sequence. Ready.", PROGRESS_NONE);
}

/* requires seq->nb_layers and seq->number to be already set */
void check_or_allocate_regparam(sequence *seq, int layer) {
	assert(layer < seq->nb_layers);
	if (!seq->regparam && seq->nb_layers > 0) {
		seq->regparam = calloc(seq->nb_layers, sizeof(regdata*));
		seq->layers = calloc(seq->nb_layers, sizeof(layer_info));
	}
	if (seq->regparam && !seq->regparam[layer] && seq->number > 0) {
		seq->regparam[layer] = calloc(seq->number, sizeof(regdata));
	}
}

/* internal sequence are a set of 1-layer images already loaded elsewhere, and
 * directly referenced as fits *.
 * This is used in LRGV composition.
 * The returned sequence does not contain any reference to files, and thus has
 * to be populated with internal_sequence_set() */
sequence *create_internal_sequence(int size) {
	int i;
	sequence *seq = calloc(1, sizeof(sequence));
	initialize_sequence(seq, TRUE);
	seq->type = SEQ_INTERNAL;
	seq->number = size;
	seq->selnum = size;
	seq->nb_layers = 1;	
	seq->internal_fits = calloc(size, sizeof(fits *));
	seq->seqname = strdup("internal sequence");
	seq->imgparam = calloc(size, sizeof(imgdata));
	for (i = 0; i < size; i++) {
		seq->imgparam[i].filenum = i;
		seq->imgparam[i].incl = 1;
		seq->imgparam[i].stats = NULL;
	}
	check_or_allocate_regparam(seq, 0);
	return seq;
}

void internal_sequence_set(sequence *seq, int index, fits *fit) {
	assert(seq);
	assert(seq->internal_fits);
	assert(index < seq->number);
	seq->internal_fits[index] = fit;
}

// find index of the fit argument in the sequence
int internal_sequence_find_index(sequence *seq, fits *fit) {
	int i;
	assert(seq);
	assert(seq->internal_fits);
	for (i = 0; i < seq->number; i++) {
		if (fit == seq->internal_fits[i])
			return i;
	}
	return -1;
}

gboolean end_crop_sequence(gpointer p) {
	struct crop_sequence_data *args = (struct crop_sequence_data *) p;
	stop_processing_thread();// can it be done here in case there is no thread?
	char *rseqname = malloc(strlen(args->prefix) + strlen(com.seq.seqname) + 5);

	sprintf(rseqname, "%s%s.seq", args->prefix, com.seq.seqname);
	check_seq(0);
	update_sequences_list(rseqname);
	set_cursor_waiting(FALSE);
	update_used_memory();
	free(args);
	free(rseqname);
	return FALSE;
}

gpointer crop_sequence(gpointer p) {
	struct crop_sequence_data *args = (struct crop_sequence_data *) p;
	int frame, ret;
	float cur_nb;

	/* then we compare to other frames */
	for (frame = 0, cur_nb = 0.f; frame < args->seq->number; frame++) {
		if (!get_thread_run())
			break;
		ret = seq_read_frame(args->seq, frame, &(wfit[0]));
		if (!ret) {
			char dest[256], filename[256];

			crop(&(wfit[0]), args->area);
			fit_sequence_get_image_filename(args->seq, frame, filename, TRUE);

			sprintf(dest, "%s%s", args->prefix, filename);
			savefits(dest, &wfit[0]);

			cur_nb += 1.f;
			set_progress_bar_data(NULL, cur_nb / args->seq->number);
		}
	}
	gdk_threads_add_idle(end_crop_sequence, args);
	return 0;
}

// check if the passed sequence is used as a color sequence. It can be a CFA
// sequence explicitly demoisaiced too, which returns true.
gboolean sequence_is_rgb(sequence *seq) {
	switch (seq->type) {
		case SEQ_REGULAR:
			return seq->nb_layers == 3;
		case SEQ_SER:
			return (seq->ser_file->color_id != SER_MONO && !com.raw_set.ser_cfa) ||
				seq->ser_file->color_id == RGB || seq->ser_file->color_id == BGR;
		default:
			return TRUE;
	}
}

/* Get statistics for an image in a sequence.
 * If it's not in the cache, it will be computed from the_image. If the_image is NULL,
 * it returns NULL in that case.
 * Do not free result.
 */
imstats* seq_get_imstats(sequence *seq, int index, fits *the_image) {
	assert(seq->imgparam);
	if (!seq->imgparam[index].stats && the_image) {
		seq->imgparam[index].stats = statistics(the_image, 0, NULL);
		seq->needs_saving = TRUE;
	}
	return seq->imgparam[index].stats;
}

