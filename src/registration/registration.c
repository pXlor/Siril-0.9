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
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <fftw3.h>
#include <assert.h>
#include <float.h>
#include <gtk/gtk.h>

#include "registration/registration.h"
#include "registration/matching/misc.h"
#include "registration/matching/match.h"
#include "algos/star_finder.h"
#include "stacking/stacking.h"
#include "core/siril.h"
#include "gui/callbacks.h"
#include "core/proto.h"
#include "algos/PSF.h"
#include "gui/PSF_list.h"
#include "algos/quality.h"
#ifdef HAVE_OPENCV
#include "opencv/opencv.h"
#endif

#define MAX_STARS_FITTED 100
#undef DEBUG

static char *tooltip_text[] = { "Image pattern alignment: This is a simple registration by translation method "
		"using cross correlation in the spatial domain. This method is fast and is used to register "
		"planetary movies. It can also be used for some deep-sky images registration. "
		"Shifts at pixel precision are saved in seq file.",
		"One star registration: This is the simplest method to register deep-sky images. "
		"Because only one star is concerned for register, images are aligned using shifting"
		"(at a fraction of pixel). No rotation or scaling are performed. "
		"Shifts at pixel precision are saved in seq file.",
		"Global star alignment: This is a more powerfull and accurate algorithm (but also slower) "
		"to perform deep-sky images. The global matching is based on triangle similarity method for automatically "
		"identify commom stars in each image."
		"A new sequence is created with the prefix of your choice (r_ by default). "
};
/* callback for the selected area event */
void _reg_selected_area_callback() {
	update_reg_interface(TRUE);
}

static struct registration_method *reg_methods[4];
static gpointer register_thread_func(gpointer p);
static gboolean end_register_idle(gpointer p);

struct registration_method *new_reg_method(char *name, registration_function f,
		selection_type s) {
	struct registration_method *reg = malloc(sizeof(struct registration_method));
	reg->name = strdup(name);
	reg->method_ptr = f;
	reg->sel = s;
	return reg;
}

void initialize_registration_methods() {
	GtkComboBoxText *regcombo;
	int i = 0, j = 0;
	GString *tip;
	gchar *ctip;

	reg_methods[i++] = new_reg_method("Image pattern alignment (planetary/deep-sky)",
			&register_shift_dft, REQUIRES_SQUARED_SELECTION);
	reg_methods[i++] = new_reg_method("One star registration (deep-sky)",
			&register_shift_fwhm, REQUIRES_ANY_SELECTION);
#ifdef HAVE_OPENCV
	reg_methods[i++] = new_reg_method("Global star alignement (deep-sky)",
			&register_star_alignment, REQUIRES_NO_SELECTION);
#endif
	//if (theli_is_available())
	//	reg_methods[i++] = new_reg_method("theli", register_theli, REQUIRES_NO_SELECTION);
	reg_methods[i] = NULL;

	tip = g_string_new ("");
	for (j = 0; j < i; j ++) {
		g_string_append(tip, tooltip_text[j]);
		if (j < i - 1)
			g_string_append(tip, "\n\n");
	}
	ctip = g_string_free (tip, FALSE);
	gtk_widget_set_tooltip_text(lookup_widget("comboboxregmethod"), ctip);
	g_free(ctip);

	/* fill comboboxregmethod */
	regcombo = GTK_COMBO_BOX_TEXT(
			gtk_builder_get_object(builder, "comboboxregmethod"));
	gtk_combo_box_text_remove_all(regcombo);
	i = 0;
	while (reg_methods[i] != NULL) {
		gtk_combo_box_text_append_text(regcombo, reg_methods[i]->name);
		siril_log_message("Added a registration method: %s\n",
				reg_methods[i]->name);
		i++;
	}
	if (i > 0) {
		gtk_combo_box_set_active(GTK_COMBO_BOX(regcombo), com.reg_settings);
	}

	/* register to the new area selected event */
	register_selection_update_callback(_reg_selected_area_callback);
}

struct registration_method *get_selected_registration_method() {
	GtkComboBox *regcombo = GTK_COMBO_BOX(
			gtk_builder_get_object(builder, "comboboxregmethod"));
	int index = gtk_combo_box_get_active(regcombo);
	return reg_methods[index];
}

/* register images: calculate shift in images to be aligned with the reference image;
 * images are not modified, only shift parameters are saved in regparam in the sequence.
 * layer is the layer on which the registration will be done, green by default (set in siril_init())
 */
int register_shift_dft(struct registration_args *args) {
	int i, x, size, sqsize;
	fftw_complex *ref, *img, *in, *out, *convol;
	fftw_plan p, q;
	int shiftx, shifty, shift;
	int ret;
	int plan;
	fftw_complex *cbuf, *ibuf, *obuf;
	float nb_frames, cur_nb;
	int ref_image;
	regdata *current_regdata;
	char tmpmsg[1024], tmpfilename[256];
	rectangle full_area;	// the area to use after getting image_part
	int j;
	double q_max = 0;
	int q_index = -1;
	fits fit;
	memset(&fit, 0, sizeof(fits));

	/* the selection needs to be squared for the DFT */
	assert(args->selection.w == args->selection.h);
	size = args->selection.w;
	full_area.x = full_area.y = 0;
	full_area.h = full_area.w = size;
	sqsize = size * size;

	if (args->process_all_frames)
		nb_frames = (float) args->seq->number;
	else
		nb_frames = (float) args->seq->selnum;

	if (!args->seq->regparam) {
		fprintf(stderr, "regparam should have been created before\n");
		return -1;
	}
	if (args->seq->regparam[args->layer]) {
		siril_log_message(
				"Recomputing already existing registration for this layers\n");
		current_regdata = args->seq->regparam[args->layer];
	} else {
		current_regdata = calloc(args->seq->number, sizeof(regdata));
		if (current_regdata == NULL) {
			siril_log_message("Error allocating registration data\n");
			return -2;
		}
	}

	/* loading reference frame */
	if (args->seq->reference_image == -1)
		ref_image = 0;
	else
		ref_image = args->seq->reference_image;

	set_progress_bar_data(
			"Register DFT: loading and processing reference frame",
			PROGRESS_NONE);
	ret = seq_read_frame_part(args->seq, args->layer, ref_image, &fit,
			&args->selection);

	q_max = current_regdata[ref_image].quality;
	q_index = ref_image;

	if (ret) {
		siril_log_message(
				"Register: could not load first image to register, aborting.\n");
		free(current_regdata);
		return ret;
	}

	ref = fftw_malloc(sizeof(fftw_complex) * sqsize);
	img = fftw_malloc(sizeof(fftw_complex) * sqsize);
	in = fftw_malloc(sizeof(fftw_complex) * sqsize);
	out = fftw_malloc(sizeof(fftw_complex) * sqsize);
	convol = fftw_malloc(sizeof(fftw_complex) * sqsize);

	if (nb_frames > 200.f)
		plan = FFTW_MEASURE;
	else
		plan = FFTW_ESTIMATE;

	p = fftw_plan_dft_2d(size, size, ref, out, FFTW_FORWARD, plan);
	q = fftw_plan_dft_2d(size, size, convol, out, FFTW_BACKWARD, plan);

	// copying image selection into the fftw data
	for (j = 0; j < sqsize; j++)
		ref[j] = (double) fit.data[j];

	// We don't need fit anymore, we can destroy it.
	current_regdata[ref_image].quality = QualityEstimate(&fit, args->layer, QUALTYPE_NORMAL);
	clearfits(&fit);
	fftw_execute_dft(p, ref, in); /* repeat as needed */
	current_regdata[ref_image].shiftx = 0;
	current_regdata[ref_image].shifty = 0;

	for (i = 0, cur_nb = 0.f; i < args->seq->number; ++i) {
		if (args->run_in_thread && !get_thread_run())
			break;
		if (i == ref_image)
			continue;
		if (!args->process_all_frames && !args->seq->imgparam[i].incl)
			continue;

		seq_get_image_filename(args->seq, i, tmpfilename);
		g_snprintf(tmpmsg, 1024, "Register: processing image %s\n", tmpfilename);
		set_progress_bar_data(tmpmsg, PROGRESS_NONE);
		if (!(ret = seq_read_frame_part(args->seq, args->layer, i, &fit,
				&args->selection))) {

			// copying image selection into the fftw data
			for (j = 0; j < sqsize; j++)
				img[j] = (double) fit.data[j];

			// We don't need fit anymore, we can destroy it.
			current_regdata[i].quality = QualityEstimate(&fit, args->layer, QUALTYPE_NORMAL);

			if (current_regdata[i].quality > q_max) {
				q_max = current_regdata[i].quality;
				q_index = i;
			}

			clearfits(&fit);
			fftw_execute_dft(p, img, out); /* repeat as needed */
			/* originally, quality is computed with the quality
			 * function working on the fft space. Now we use quality
			 * instead.
			 */
			cbuf = convol;
			ibuf = in;
			obuf = out;
			for (x = 0; x < sqsize; x++) {
				*cbuf++ = *ibuf++ * conj(*obuf++);
			}

			fftw_execute_dft(q, convol, ref); /* repeat as needed */
			shift = 0;
			for (x = 1; x < sqsize; ++x) {
				if (creal(ref[x]) > creal(ref[shift])) {
					shift = x;
					// break or get last value?
				}
			}
			shifty = shift / size;
			shiftx = shift % size;
			if (shifty > size / 2) {
				shifty -= size;
			}
			if (shiftx > size / 2) {
				shiftx -= size;
			}

			current_regdata[i].shiftx = shiftx;
			current_regdata[i].shifty = shifty;

			/* shiftx and shifty are the x and y values for translation that
			 * would make this image aligned with the reference image.
			 * WARNING: the y value is counted backwards, since the FITS is
			 * stored down from up.
			 */
			fprintf(stderr, "reg: file %d, shiftx=%d shifty=%d quality=%g\n",
					args->seq->imgparam[i].filenum, current_regdata[i].shiftx,
					current_regdata[i].shifty, current_regdata[i].quality);
			cur_nb += 1.f;
			set_progress_bar_data(NULL, cur_nb / nb_frames);
		} else {
			//report_fits_error(ret, error_buffer);
			if (current_regdata == args->seq->regparam[args->layer])
				args->seq->regparam[args->layer] = NULL;
			free(current_regdata);
			return ret;
		}
	}

	fftw_destroy_plan(p);
	fftw_destroy_plan(q);
	fftw_free(in);
	fftw_free(out);
	fftw_free(ref);
	fftw_free(img);
	fftw_free(convol);
	args->seq->regparam[args->layer] = current_regdata;
	update_used_memory();
	siril_log_message("Registration finished.\n");
	siril_log_color_message("Best frame: #%d with quality=%g.\n", "bold",
			q_index, q_max);
	return 0;
}

/* register images: calculate shift in images to be aligned with the reference image;
 * images are not modified, only shift parameters are saved in regparam in the sequence.
 * layer is the layer on which the registration will be done, green by default (set in siril_init())
 */
int register_shift_fwhm(struct registration_args *args) {
	int frame, ref_image;
	float nb_frames, cur_nb;
	double reference_xpos, reference_ypos;
	double fwhm_min = DBL_MAX;
	int fwhm_index = -1;
	regdata *current_regdata;
	/* First and longest step: get the minimization data on one star for all
	 * images to register, which provides FWHM but also star coordinates */
	// TODO: detect that it was already computed, and don't do it again
	// -> should be done at a higher level and passed in the args
	if (do_fwhm_sequence_processing(args->seq, args->layer))	// stores in regparam
		return 1;

	current_regdata = args->seq->regparam[args->layer];

	if (args->process_all_frames)
		nb_frames = (float) args->seq->number;
	else
		nb_frames = (float) args->seq->selnum;

	/* loading reference frame */
	if (args->seq->reference_image == -1)
		ref_image = 0;
	else
		ref_image = args->seq->reference_image;
	if (!current_regdata[ref_image].fwhm_data) {
		siril_log_message(
				"Registration PSF: failed to compute PSF for reference frame at least\n");
		if (current_regdata != args->seq->regparam[args->layer])
			free(current_regdata);
		return -1;
	}
	reference_xpos = current_regdata[ref_image].fwhm_data->xpos;
	reference_ypos = current_regdata[ref_image].fwhm_data->ypos;

	fwhm_min = current_regdata[ref_image].fwhm_data->fwhmx;

	fwhm_index = ref_image;

	/* Second step: align image by aligning star coordinates together */
	for (frame = 0, cur_nb = 0.f; frame < args->seq->number; frame++) {
		double tmp;
		if (args->run_in_thread && !get_thread_run())
			break;
		if (!args->process_all_frames && !args->seq->imgparam[frame].incl)
			continue;
		if (frame == ref_image || !current_regdata[frame].fwhm_data) {
			current_regdata[frame].shiftx = 0;
			current_regdata[frame].shifty = 0;
			continue;
		}
		if (current_regdata[frame].fwhm < fwhm_min
				&& current_regdata[frame].fwhm > 0.0) {
			fwhm_min = current_regdata[frame].fwhm;
			fwhm_index = frame;
		}
		tmp = reference_xpos - current_regdata[frame].fwhm_data->xpos;
		current_regdata[frame].shiftx = round_to_int(tmp);
		tmp = current_regdata[frame].fwhm_data->ypos - reference_ypos;
		current_regdata[frame].shifty = round_to_int(tmp);

		/* shiftx and shifty are the x and y values for translation that
		 * would make this image aligned with the reference image.
		 * WARNING: the y value is counted backwards, since the FITS is
		 * stored down from up.
		 */
		fprintf(stderr, "reg: file %d, shiftx=%d shifty=%d\n",
				args->seq->imgparam[frame].filenum,
				current_regdata[frame].shiftx, current_regdata[frame].shifty);
		cur_nb += 1.f;
		set_progress_bar_data(NULL, cur_nb / nb_frames);
	}

	args->seq->regparam[args->layer] = current_regdata;
	update_used_memory();
	siril_log_message("Registration finished.\n");
	siril_log_color_message("Best frame: #%d with fwhm=%.3g.\n", "bold",
			fwhm_index, fwhm_min);
	return 0;
}

#ifdef HAVE_OPENCV
void _print_result(TRANS *trans, float FWHMx, float FWHMy) {
	double rotation, scale;
	point shift;

	switch (trans->order) {
	case 1:
		rotation = atan2(trans->c, trans->b);
		shift.x = trans->a;
		shift.y = -trans->d;
		scale = sqrt(trans->b * trans->b + trans->c * trans->c);
		siril_log_color_message("Matching stars: done\n", "green");
		siril_log_message("%d pair matches.\n", trans->nr);
		siril_log_message("scale:\t%*.3f\n", 9, scale);
		siril_log_message("rotation:\t%*.2f deg\n", 8, rotation * 180 / M_PI);
		siril_log_message("dx:\t\t%*.2f px\n", 8, shift.x);
		siril_log_message("dy:\t\t%*.2f px\n", 8, shift.y);
		siril_log_message("FWHMx:\t%*.2f px\n", 8, FWHMx);
		siril_log_message("FWHMy:\t%*.2f px\n", 8, FWHMy);
		break;
	default:
		siril_log_color_message("Not handled yet\n", "red");
	}
}

int register_star_alignment(struct registration_args *args) {
	int frame, ref_image, ret, i;
	int fitted_stars;
	float nb_frames, cur_nb;
	float FWHMx, FWHMy;
	fitted_PSF **stars;
	TRANS trans;
	regdata *current_regdata;
	starFinder sf;
	fits fit;
	memset(&fit, 0, sizeof(fits));

	memset(&sf, 0, sizeof(starFinder));

	if (!args->seq->regparam) {
		fprintf(stderr, "regparam should have been created before\n");
		return -1;
	}
	if (args->seq->regparam[args->layer]) {
		current_regdata = args->seq->regparam[args->layer];
		/* we reset all values as we built another sequence */
		memset(current_regdata, 0, args->seq->number * sizeof(regdata));
	} else {
		current_regdata = calloc(args->seq->number, sizeof(regdata));
		if (current_regdata == NULL) {
			siril_log_message("Error allocating registration data\n");
			return -2;
		}
	}

	if (args->process_all_frames)
		nb_frames = (float) args->seq->number;
	else
		nb_frames = (float) args->seq->selnum;

	/* loading reference frame */
	if (args->seq->reference_image == -1)
		ref_image = 0;
	else
		ref_image = args->seq->reference_image;

	/* first we're looking for stars in reference image */
	ret = seq_read_frame(args->seq, ref_image, &fit);
	if (ret) {
		siril_log_message("Could not load reference image\n");
		if (current_regdata == args->seq->regparam[args->layer])
			args->seq->regparam[args->layer] = NULL;
		free(current_regdata);
		return 1;
	}
	siril_log_color_message("Reference Image:\n", "green");
	com.stars = peaker(&fit, args->layer, &sf);
	if (sf.nb_stars < 4) {
		siril_log_message(
				"There are not enough stars in reference image to perform alignment\n");
		if (current_regdata == args->seq->regparam[args->layer])
			args->seq->regparam[args->layer] = NULL;
		free(current_regdata);
		return 1;
	}
	redraw(com.cvport, REMAP_NONE); // draw stars
#ifdef DEBUG
	printf("REFERENCE IMAGE\n");
	for (i = 0; i < MAX_STARS_FITTED; i++) {
		printf("%.3lf\t%.3lf\t%.3lf\n",
				com.stars[i]->xpos, com.stars[i]->ypos, com.stars[i]->mag);
	}
#endif
	fitted_stars = (sf.nb_stars > MAX_STARS_FITTED) ? MAX_STARS_FITTED : sf.nb_stars;
	FWHM_average(com.stars, &FWHMx, &FWHMy, fitted_stars);
	siril_log_message("FWHMx:\t%*.2f px\n", 8, FWHMx);
	siril_log_message("FWHMy:\t%*.2f px\n", 8, FWHMy);
	current_regdata[ref_image].fwhm = FWHMx;

	/* then we compare to other frames */
	args->seq->new_total = args->seq->number;
	for (frame = 0, cur_nb = 0.f; frame < args->seq->number; frame++) {
		if (args->run_in_thread && !get_thread_run())
			break;
		if (!args->process_all_frames && !args->seq->imgparam[frame].incl)
			continue;

		ret = seq_read_frame(args->seq, frame, &fit);
		if (!ret) {
			char dest[256], filename[256];
			int nbpoints;

			if (frame != ref_image) {
				stars = peaker(&fit, args->layer, &sf);
				if (sf.nb_stars < 4) {
					siril_log_message("Not enough stars. Image %d skipped\n",
							frame);
					args->seq->new_total--;
					continue;
				}

#ifdef DEBUG
				printf("IMAGE %d\n", frame);
				for (i = 0; i < MAX_STARS_FITTED; i++) {
					printf("%.3lf\t%.3lf\t%.3lf\n", stars[i]->xpos,
							stars[i]->ypos, stars[i]->mag);
				}
				printf("\n\n", frame);
#endif

				nbpoints = (sf.nb_stars < fitted_stars) ?
								sf.nb_stars : fitted_stars;

				if (star_match(stars, com.stars, nbpoints, &trans)) {
					siril_log_color_message(
							"Cannot perform star matching. Image %d skipped\n",
							"red", frame);
					args->seq->new_total--;
					i = 0;
					while (i < MAX_STARS && stars[i])
						free(stars[i++]);
					free(stars);
					continue;
				}

				FWHM_average(stars, &FWHMx, &FWHMy, nbpoints);
				_print_result(&trans, FWHMx, FWHMy);
				current_regdata[frame].fwhm = FWHMx;

				transforme_image(&fit, trans, 0);

				i = 0;
				while (i < MAX_STARS && stars[i])
					free(stars[i++]);
				free(stars);
			}
			fit_sequence_get_image_filename(args->seq, frame, filename, TRUE);

			snprintf(dest, 255, "%s%s", args->text, filename);
			savefits(dest, &fit);

			cur_nb += 1.f;
			set_progress_bar_data(NULL, cur_nb / nb_frames);
		}
	}
	args->seq->regparam[args->layer] = current_regdata;
	update_used_memory();
	siril_log_message("Registration finished.\n");

	return 0;
}
#endif

void on_comboboxregmethod_changed(GtkComboBox *box, gpointer user_data) {
	int reg = gtk_combo_box_get_active(box);

	com.reg_settings = reg;
	update_reg_interface(TRUE);
	writeinitfile();
}

int get_registration_layer(sequence *seq) {
	int reglayer;
	gboolean has_changed = FALSE;
	GtkComboBox *registbox = GTK_COMBO_BOX(lookup_widget("comboboxreglayer"));

	if (!sequence_is_loaded())
		return -1;
	reglayer = gtk_combo_box_get_active(registbox);

	if (seq->nb_layers == 3) {
		if (!seq->regparam[reglayer]) {
			if (seq->regparam[GLAYER]) {
				reglayer = GLAYER;
				has_changed = TRUE;
			}
			else if (seq->regparam[RLAYER]) {
				reglayer = RLAYER;
				has_changed = TRUE;
			}
			else if (seq->regparam[BLAYER]) {
				reglayer = BLAYER;
				has_changed = TRUE;
			}
		}
		if (has_changed) {
			gtk_combo_box_set_active(registbox, reglayer);
		}
	}
	return reglayer;
}

/* Selects the "register all" or "register selected" according to the number of
 * selected images, if argument is false.
 * Verifies that enough images are selected and an area is selected.
 */
void update_reg_interface(gboolean dont_change_reg_radio) {
	static GtkWidget *go_register = NULL, *w_entry = NULL;
	static GtkLabel *labelreginfo = NULL;
	static GtkToggleButton *reg_all = NULL, *reg_sel = NULL;
	int nb_images_reg; /* the number of images to register */
	struct registration_method *method;

	if (!go_register) {
		go_register = lookup_widget("goregister_button");
		reg_all = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "regallbutton"));
		reg_sel = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "regselbutton"));
		labelreginfo = GTK_LABEL(
				gtk_builder_get_object(builder, "labelregisterinfo"));
		w_entry = lookup_widget("regseqname_entry");
	}

	if (!dont_change_reg_radio) {
		if (com.seq.selnum < com.seq.number)
			gtk_toggle_button_set_active(reg_sel, TRUE);
		else
			gtk_toggle_button_set_active(reg_all, TRUE);
	}

	/* getting the selected registration method */
	method = get_selected_registration_method();

	if (gtk_toggle_button_get_active(reg_all))
		nb_images_reg = com.seq.number;
	else
		nb_images_reg = com.seq.selnum;
	if (method && ((nb_images_reg > 1 && com.selection.w > 0 && com.selection.h > 0)
			|| (nb_images_reg > 1 && method->sel == REQUIRES_NO_SELECTION))) {
		gtk_widget_set_sensitive(go_register, TRUE);
		gtk_label_set_text(labelreginfo, "");
#ifdef HAVE_OPENCV
		gtk_widget_set_sensitive(w_entry, method->method_ptr == &register_star_alignment);
#endif
	} else {
		gtk_widget_set_sensitive(go_register, FALSE);
		gtk_widget_set_sensitive(w_entry, FALSE);
		if (nb_images_reg <= 1 && com.selection.w <= 0 && com.selection.h <= 0)
			if (!sequence_is_loaded())
				gtk_label_set_text(labelreginfo, "Load a sequence first.");
			else
				gtk_label_set_text(labelreginfo,
					"Select an area in image first, and select images in the sequence.");
		else if (nb_images_reg <= 1)
			gtk_label_set_text(labelreginfo, "Select images in the sequence.");
		else
			gtk_label_set_text(labelreginfo, "Select an area in image first.");
	}
}

/* try to maximize the area within the image size (based on gfit) */
void compute_squared_selection(rectangle *area) {
	//fprintf(stdout, "function entry: %d,%d,\t%dx%d\n", area->x, area->y, area->w, area->h);
	if (area->x >= 0 && area->x + area->w <= gfit.rx && area->y >= 0
			&& area->y + area->h <= gfit.ry)
		return;

	if (area->x < 0) {
		area->x++;
		if (area->x + area->w > gfit.rx) {
			/* reduce area */
			area->w -= 2;
			area->h -= 2;
			area->y++;
		}
	} else if (area->x + area->w > gfit.rx) {
		area->x--;
		if (area->x < 0) {
			/* reduce area */
			area->x++;
			area->w -= 2;
			area->h -= 2;
			area->y++;
		}
	}

	if (area->y < 0) {
		area->y++;
		if (area->y + area->h > gfit.ry) {
			/* reduce area */
			area->h -= 2;
			area->w -= 2;
			area->x++;
		}
	} else if (area->y + area->h > gfit.ry) {
		area->y--;
		if (area->y < 0) {
			/* reduce area */
			area->x++;
			area->w -= 2;
			area->h -= 2;
			area->y++;
		}
	}

	return compute_squared_selection(area);
}

void get_the_registration_area(struct registration_args *reg_args,
		struct registration_method *method) {
	int max;
	switch (method->sel) {
	case REQUIRES_NO_SELECTION:
		break;
	case REQUIRES_ANY_SELECTION:
		memcpy(&reg_args->selection, &com.selection, sizeof(rectangle));
		break;
	case REQUIRES_SQUARED_SELECTION:
		/* Passed arguments are X,Y of the center of the square and the size of
		 * the square. */
		if (com.selection.w > com.selection.h)
			max = com.selection.w;
		else
			max = com.selection.h;

		reg_args->selection.x = com.selection.x + com.selection.w / 2 - max / 2;
		reg_args->selection.w = max;
		reg_args->selection.y = com.selection.y + com.selection.h / 2 - max / 2;
		reg_args->selection.h = max;
		compute_squared_selection(&reg_args->selection);

		/* save it back to com.selection do display it properly */
		memcpy(&com.selection, &reg_args->selection, sizeof(rectangle));
		fprintf(stdout, "final area: %d,%d,\t%dx%d\n", reg_args->selection.x,
				reg_args->selection.y, reg_args->selection.w,
				reg_args->selection.h);
		redraw(com.cvport, REMAP_NONE);
		break;
	}
}

/* callback for 'Go register' button, GTK thread */
void on_seqregister_button_clicked(GtkButton *button, gpointer user_data) {
	struct registration_args *reg_args;
	struct registration_method *method;
	char *msg;
	GtkToggleButton *regall;
	GtkComboBox *cbbt_layers;

	if (get_thread_run()) {
		siril_log_message(
				"Another task is already in progress, ignoring new request.\n");
		return;
	}

	/* getting the selected registration method */
	method = get_selected_registration_method();

	if (com.selection.w <= 0 && com.selection.h <= 0
			&& method->sel != REQUIRES_NO_SELECTION) {
		msg = siril_log_message(
						"All prerequisites are not filled for registration. Select a rectangle first.\n");
		show_dialog(msg, "Warning", "gtk-dialog-warning");
		return;
	}
	// TODO: check for reentrance

	reg_args = malloc(sizeof(struct registration_args));

	control_window_switch_to_tab(OUTPUT_LOGS);

	/* filling the arguments for registration */
	reg_args->seq = &com.seq;
	regall = GTK_TOGGLE_BUTTON(lookup_widget("regallbutton"));
	reg_args->process_all_frames = gtk_toggle_button_get_active(regall);
	/* getting the selected registration layer from the combo box. The value is the index
	 * of the selected line, and they are in the same order than layers so there should be
	 * an exact matching between the two */
	cbbt_layers = GTK_COMBO_BOX(
			gtk_builder_get_object(builder, "comboboxreglayer"));
	reg_args->layer = gtk_combo_box_get_active(cbbt_layers);
	get_the_registration_area(reg_args, method);
	reg_args->func = method->method_ptr;
	reg_args->run_in_thread = TRUE;
	reg_args->entry = GTK_ENTRY(
			gtk_builder_get_object(builder, "regseqname_entry"));
	reg_args->text = gtk_entry_get_text(reg_args->entry);

	msg = siril_log_color_message("Registration: processing using method: %s\n",
			"red", method->name);
	msg[strlen(msg) - 1] = '\0';
	gettimeofday(&(reg_args->t_start), NULL);
	set_cursor_waiting(TRUE);
	set_progress_bar_data(msg, PROGRESS_RESET);

	start_in_new_thread(register_thread_func, reg_args);
}

// worker thread
static gpointer register_thread_func(gpointer p) {
	struct registration_args *args = (struct registration_args *) p;
	args->retval = args->func(args);
	gdk_threads_add_idle(end_register_idle, args);
	return GINT_TO_POINTER(args->retval);	// not used anyway
}

// end of registration, GTK thread
static gboolean end_register_idle(gpointer p) {
	struct timeval t_end;
	struct registration_args *args = (struct registration_args *) p;
	stop_processing_thread();
	if (!args->retval) {
		writeseqfile(args->seq);
		fill_sequence_list(args->seq, com.cvport);
		set_layers_for_registration();	// update display of available reg data

		/* Load new sequence. Only star alignment method uses new sequence. */
#ifdef HAVE_OPENCV
		if (args->func == &register_star_alignment) {
			int frame, new_frame;
			regdata *new_data;
			imgdata *new_image;
			char *rseqname = malloc(gtk_entry_get_text_length(args->entry)
							+ strlen(com.seq.seqname) + 5);

			sprintf(rseqname, "%s%s.seq", args->text, com.seq.seqname);
			unlink(rseqname);
			check_seq(0);
			if (args->seq->seqname)
				free(args->seq->seqname);
			char *newname = remove_ext_from_filename(rseqname);
			args->seq->seqname = strdup(newname);

			/* If images have not been registred we have to reorganize data and images */
			new_data = calloc(args->seq->new_total, sizeof(regdata));
			new_image = calloc(args->seq->new_total, sizeof(imgdata));
			for (frame = 0, new_frame = 0; frame < args->seq->number; frame++) {
				if (!args->process_all_frames
						&& !args->seq->imgparam[frame].incl)
					continue;
				if (args->seq->regparam[args->layer][frame].fwhm > 0) {
					new_data[new_frame] = args->seq->regparam[args->layer][frame];
					new_image[new_frame] = args->seq->imgparam[frame];
					new_frame++;
				}
			}
			if (args->seq->regparam[args->layer])
				free(args->seq->regparam[args->layer]);
			if (args->seq->imgparam)
				free(args->seq->imgparam);
			args->seq->regparam[args->layer] = new_data;
			args->seq->imgparam = new_image;
			args->seq->number = args->seq->new_total;
			args->seq->selnum = args->seq->new_total;

			/* We write the new sequence */
			writeseqfile(args->seq);
			update_sequences_list(rseqname);
			free(newname);
			free(rseqname);
			clear_stars_list();
		}
#endif
	}
	set_progress_bar_data("Registration complete", PROGRESS_DONE);
	update_stack_interface();
	set_cursor_waiting(FALSE);
	gettimeofday(&t_end, NULL);
	show_time(args->t_start, t_end);
	free(args);
	return FALSE;
}
