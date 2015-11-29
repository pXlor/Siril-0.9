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

#include <gsl/gsl_histogram.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>
#include "core/siril.h"
#include "core/proto.h"
#include "io/single_image.h"
#include "gui/histogram.h"
#include "gui/callbacks.h"	// for lookup_widget()
#include "core/undo.h"

/* The gsl_histogram, documented here:
 * http://linux.math.tifr.res.in/manuals/html/gsl-ref-html/gsl-ref_21.html
 * is able to group values into bins, it does not need to handle all values
 * separately. That is useful for display purposes, but not currently used.
 */

// margin between axis and drawing area border
//#define AXIS_MARGIN 10
// colors of layers histograms		R	G	B	RGB
static double histo_color_r[] = { 1.0, 0.0, 0.0, 0.0 };
static double histo_color_g[] = { 0.0, 1.0, 0.0, 0.0 };
static double histo_color_b[] = { 0.0, 0.0, 1.0, 0.0 };
static double graph_height = 0.;
static uint64_t clipped[] = { 0, 0 };

static GtkToggleButton *toggles[MAXVPORT] = { NULL };

static void _init_clipped_pixels() {
	clipped[0] = 0;
	clipped[1] = 0;
}

static void _initialize_clip_text() {
	static GtkEntry *clip_high = NULL, *clip_low = NULL;

	if (clip_high == NULL) {
		clip_high = GTK_ENTRY(
				gtk_builder_get_object(builder, "clip_highlights"));
		clip_low = GTK_ENTRY(gtk_builder_get_object(builder, "clip_shadows"));
	}
	gtk_entry_set_text(clip_low, "0.000%");
	gtk_entry_set_text(clip_high, "0.000%");
}

static void _update_clipped_pixels(int data) {
	static GtkEntry *clip_high = NULL, *clip_low = NULL;
	double tmp;
	char buffer[16];

	if (clip_high == NULL) {
		clip_high = GTK_ENTRY(lookup_widget("clip_highlights"));
		clip_low = GTK_ENTRY(lookup_widget("clip_shadows"));
	}
	tmp = (double) clipped[1] * 100.0 / data;
	g_snprintf(buffer, sizeof(buffer), "%.3lf%%", tmp);
	gtk_entry_set_text(clip_high, buffer);
	tmp = (double) clipped[0] * 100.0 / data;
	g_snprintf(buffer, sizeof(buffer), "%.3lf%%", tmp);
	gtk_entry_set_text(clip_low, buffer);

}

gsl_histogram* computeHisto_Selection(fits* fit, int layer,
		rectangle *selection) {
	WORD *from;
	size_t size = (size_t) get_normalized_value(fit);
	gsl_histogram* histo = gsl_histogram_alloc(size + 1);
	int stridefrom, i, j;

	gsl_histogram_set_ranges_uniform(histo, 0, size);
	from = fit->pdata[layer] + (fit->ry - selection->y - selection->h) * fit->rx
			+ selection->x;
	stridefrom = fit->rx - selection->w;
	for (i = 0; i < selection->h; i++) {
		for (j = 0; j < selection->w; j++) {
			gsl_histogram_increment(histo, (double) *from);
			from++;
		}
		from += stridefrom;
	}
	return histo;
}

gsl_histogram* computeHisto(fits* fit, int layer) {
	unsigned int i, ndata;
	WORD *buf;
	if (layer >= 3)
		return NULL;
	size_t size = (size_t) get_normalized_value(fit) + 1;
	gsl_histogram* histo = gsl_histogram_alloc(size);
	gsl_histogram_set_ranges_uniform(histo, 0, size - 1);

	buf = fit->pdata[layer];
	ndata = fit->rx * fit->ry;
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)
	for (i = 0; i < ndata; i++) {
		gsl_histogram_increment(histo, (double) buf[i]);
	}
	return histo;
}

//compute histogram for all pixel values below backgroud value
gsl_histogram* histo_bg(fits* fit, int layer, double bg) {
	unsigned int i, ndata;
	WORD *buf;
	if (layer >= 3)
		return NULL;
	size_t size = (size_t) bg;
	gsl_histogram* histo = gsl_histogram_alloc(size + 1);
	gsl_histogram_set_ranges_uniform(histo, 0, size);

	buf = fit->pdata[layer];
	ndata = fit->rx * fit->ry;
	for (i = 0; i < ndata; i++) {
		if (buf[i] <= (WORD) bg)
			gsl_histogram_increment(histo, (double) buf[i]);
	}
	return histo;
}

void compute_histo_for_gfit(int force) {
	int nb_layers = 3;
	int i;
	if (gfit.naxis == 2)
		nb_layers = 1;
	for (i = 0; i < nb_layers; i++) {
		if (force || !com.layers_hist[i])
			set_histogram(computeHisto(&gfit, i), i);
	}
	set_histo_toggles_names();
}

void update_gfit_histogram_if_needed() {
	static GtkWidget *drawarea = NULL, *selarea = NULL;
	if (is_histogram_visible())
		compute_histo_for_gfit(1);
	if (!drawarea) {
		drawarea = lookup_widget("drawingarea_histograms");
		selarea = lookup_widget("drawingarea_histograms_selection");
	}
	gtk_widget_queue_draw(drawarea);
	gtk_widget_queue_draw(selarea);
}

void _histo_on_selection_changed() {
	static GtkWidget *selarea = NULL;
	if (!selarea)
		selarea = lookup_widget("drawingarea_histograms_selection");
	gtk_widget_queue_draw(selarea);
}

void set_histogram(gsl_histogram *histo, int layer) {
	assert(layer >= 0 && layer < MAXVPORT);
	if (com.layers_hist[layer])
		gsl_histogram_free(com.layers_hist[layer]);
	com.layers_hist[layer] = histo;
}

void clear_histograms() {
	int i;
	for (i = 0; i < MAXVPORT; i++)
		set_histogram(NULL, i);
}

void init_toggles() {
	if (!toggles[0]) {
		toggles[0] = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "check_histo_r"));
		toggles[1] = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "check_histo_g"));
		toggles[2] = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "check_histo_b"));
		toggles[3] = NULL; //GTK_TOGGLE_BUTTON(gtk_builder_get_object(builder, "check_histo_rgb"));
	}
}

// sets the channel names of the toggle buttons in the histogram window, based on
// the number of layers of gfit
void set_histo_toggles_names() {
	init_toggles();
	if (gfit.naxis == 2) {
		const char* test = gtk_button_get_label(GTK_BUTTON(toggles[0]));
		if (strcmp(test, "gray"))
			gtk_button_set_label(GTK_BUTTON(toggles[0]), "gray");
		gtk_toggle_button_set_active(toggles[0], TRUE);
		gtk_widget_set_visible(GTK_WIDGET(toggles[1]), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(toggles[2]), FALSE);
		/* visible has no effect in GTK+ 3.12, trying sensitive too 
		 * Yes it does. The solution is to call the window (widget)
		 * with gtk_widget_show and not gtk_widget_show_all */
		gtk_widget_set_sensitive(GTK_WIDGET(toggles[1]), FALSE);
		gtk_widget_set_sensitive(GTK_WIDGET(toggles[2]), FALSE);
		if (toggles[3])
			gtk_widget_set_visible(GTK_WIDGET(toggles[3]), FALSE);
	} else {
		const char* test = gtk_button_get_label(GTK_BUTTON(toggles[0]));
		if (strcmp(test, "red"))
			gtk_button_set_label(GTK_BUTTON(toggles[0]), "red");
		gtk_toggle_button_set_active(toggles[0], TRUE);
		gtk_toggle_button_set_active(toggles[1], TRUE);
		gtk_toggle_button_set_active(toggles[2], TRUE);
		gtk_widget_set_sensitive(GTK_WIDGET(toggles[1]), TRUE);
		gtk_widget_set_sensitive(GTK_WIDGET(toggles[2]), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(toggles[1]), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(toggles[2]), TRUE);
		if (toggles[3]) {
			gtk_widget_set_visible(GTK_WIDGET(toggles[3]), TRUE);
			gtk_toggle_button_set_active(toggles[3], TRUE);
		}
	}
}

gboolean redraw_histo(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int i, width, height;

	fprintf(stdout, "histogram redraw\n");
	init_toggles();
	width = gtk_widget_get_allocated_width(widget);
	height = gtk_widget_get_allocated_height(widget);
	if (height == 1)
		return FALSE;
	erase_histo_display(cr, width, height);
	graph_height = 0.0;
	for (i = 0; i < MAXVPORT; i++) {
		if (com.layers_hist[i]
				&& (!toggles[i] || gtk_toggle_button_get_active(toggles[i])))
			display_histo(com.layers_hist[i], cr, i, width, height);
	}
	return FALSE;
}

gboolean redraw_histo_selection(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int i, width, height, nb_layers;
	double tmp_height = graph_height;

	width = gtk_widget_get_allocated_width(widget);
	height = gtk_widget_get_allocated_height(widget);
	if (height == 1)
		return FALSE;
	erase_histo_display(cr, width, height);

	if (!com.selection.w || !com.selection.h)
		return TRUE;
	if (sequence_is_loaded())
		nb_layers = com.seq.nb_layers;
	else
		nb_layers = com.uniq->nb_layers;

	for (i = 0; i < nb_layers; i++) {
		gsl_histogram *histo = computeHisto_Selection(&gfit, i, &com.selection);
		graph_height = 0.0;
		if (histo)
			display_histo(histo, cr, i, width, height);
	}
	graph_height = tmp_height;
	return FALSE;
}

void on_histo_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	static GtkWidget *drawarea = NULL, *selarea = NULL;
	if (!drawarea) {
		drawarea = lookup_widget("drawingarea_histograms");
		selarea = lookup_widget("drawingarea_histograms_selection");
	}
	gtk_widget_queue_draw(drawarea);
	gtk_widget_queue_draw(selarea);
}

// erase image and redraw the background color and grid
void erase_histo_display(cairo_t *cr, int width, int height) {
	double dash_format[] = { 8.0, 4.0 };
	// clear all with background color
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_rectangle(cr, 0, 0, width, height);
	cairo_fill(cr);
	// draw axis
	/*cairo_set_line_width(cr, 1.0);
	 cairo_move_to(cr, AXIS_MARGIN, AXIS_MARGIN);
	 cairo_line_to(cr, AXIS_MARGIN, height-AXIS_MARGIN);
	 cairo_line_to(cr, width-AXIS_MARGIN, height-AXIS_MARGIN);*/
	// draw grid 
	cairo_set_line_width(cr, 1.0);
	cairo_set_source_rgb(cr, 0.4, 0.4, 0.4);
	// quarters in solid, eights in dashed line
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_move_to(cr, width * 0.25, 0);
	cairo_line_to(cr, width * 0.25, height);
	cairo_move_to(cr, width * 0.5, 0);
	cairo_line_to(cr, width * 0.5, height);
	cairo_move_to(cr, width * 0.75, 0);
	cairo_line_to(cr, width * 0.75, height);

	cairo_move_to(cr, 0, height * 0.25);
	cairo_line_to(cr, width, height * 0.25);
	cairo_move_to(cr, 0, height * 0.5);
	cairo_line_to(cr, width, height * 0.5);
	cairo_move_to(cr, 0, height * 0.75);
	cairo_line_to(cr, width, height * 0.75);

	cairo_stroke(cr);

	cairo_set_dash(cr, dash_format, 2, 0);
	cairo_move_to(cr, width * 0.125, 0);
	cairo_line_to(cr, width * 0.125, height);
	cairo_move_to(cr, width * 0.375, 0);
	cairo_line_to(cr, width * 0.375, height);
	cairo_move_to(cr, width * 0.625, 0);
	cairo_line_to(cr, width * 0.625, height);
	cairo_move_to(cr, width * 0.875, 0);
	cairo_line_to(cr, width * 0.875, height);

	cairo_move_to(cr, 0, height * 0.125);
	cairo_line_to(cr, width, height * 0.125);
	cairo_move_to(cr, 0, height * 0.375);
	cairo_line_to(cr, width, height * 0.375);
	cairo_move_to(cr, 0, height * 0.625);
	cairo_line_to(cr, width, height * 0.625);
	cairo_move_to(cr, 0, height * 0.875);
	cairo_line_to(cr, width, height * 0.875);
	cairo_stroke(cr);
}

void display_histo(gsl_histogram *histo, cairo_t *cr, int layer, int width,
		int height) {
	int current_bin;
	size_t norm = gsl_histogram_bins(histo) - 1;

	float vals_per_px = (float) (norm) / (float) width;	// size of a bin
	size_t i, nb_orig_bins = gsl_histogram_bins(histo);

	// We need to store the binned histogram in order to find the binned maximum
	static double *displayed_values = NULL;
	static int nb_bins_allocated = 0;
	/* we create a bin for each pixel in the displayed width.
	 * nb_bins_allocated is thus equal to the width of the image */
	if (nb_bins_allocated != width) {
		double *tmp;
		nb_bins_allocated = width;
		tmp = realloc(displayed_values, nb_bins_allocated * sizeof(double));
		if (!tmp) {
			if (displayed_values)
				free(displayed_values);
			fprintf(stderr, "Failed to reallocate histogram bins\n");
			return;
		}
		displayed_values = tmp;
		memset(displayed_values, 0, nb_bins_allocated);
	}
	assert(displayed_values);

	if (gfit.naxis == 2)
		cairo_set_source_rgb(cr, 255.0, 255.0, 255.0);
	else
		cairo_set_source_rgb(cr, histo_color_r[layer], histo_color_g[layer],
				histo_color_b[layer]);
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_set_line_width(cr, 1.5);

	// first loop builds the bins and finds the maximum
	i = 0;
	current_bin = 0;
	do {
		double bin_val = 0.0;
		while (i < nb_orig_bins
				&& (float) i / vals_per_px <= (float) current_bin + 0.5f) {
			bin_val += gsl_histogram_get(histo, i);
			i++;
		}
		displayed_values[current_bin] = bin_val;
		if (bin_val > graph_height)	// check for maximum
			graph_height = bin_val;
		current_bin++;
	} while (i < nb_orig_bins && current_bin < nb_bins_allocated);

	//graph_height *= 1.05;	// margin, but grid is not aligned anymore
	//cairo_move_to(cr, 0, height);	// first line_to will act as move_to
	for (i = 0; i < nb_bins_allocated; i++) {
		double bin_height = height
				- height * displayed_values[i] / graph_height;
		cairo_line_to(cr, i, bin_height);
	}
	cairo_stroke(cr);
}

void on_histogram_window_show(GtkWidget *object, gpointer user_data) {
	register_selection_update_callback(_histo_on_selection_changed);
	_initialize_clip_text();
}

void on_histogram_window_hide(GtkWidget *object, gpointer user_data) {
	unregister_selection_update_callback(_histo_on_selection_changed);
}

void on_button_histo_close_clicked() {
	graph_height = 0.;
	gtk_widget_hide(lookup_widget("histogram_window"));
	reset_curors_and_values();
}

void reset_curors_and_values() {
	gtk_range_set_value(
			GTK_RANGE(gtk_builder_get_object(builder, "scale_midtones")), 0.5);
	gtk_range_set_value(
			GTK_RANGE(gtk_builder_get_object(builder, "scale_shadows")), 0.0);
	gtk_range_set_value(
			GTK_RANGE(gtk_builder_get_object(builder, "scale_highlights")),
			1.0);
	_init_clipped_pixels();
	_initialize_clip_text();
	update_gfit_histogram_if_needed();
}

void on_button_histo_reset_clicked() {
	set_cursor_waiting(TRUE);
	reset_curors_and_values();
	set_cursor_waiting(FALSE);
}

void apply_mtf_to_fits(fits *fit) {
	double m, lo, hi, pente;
	int i, chan, nb_chan, ndata;
	WORD *buf[3] =
			{ fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	static GtkRange *scale_transfert_function[3] = { NULL, NULL, NULL };
	WORD norm = get_normalized_value(fit);

	if (scale_transfert_function[1] == NULL) {
		scale_transfert_function[0] = GTK_RANGE(lookup_widget("scale_shadows"));
		scale_transfert_function[1] = GTK_RANGE(
				lookup_widget("scale_midtones"));
		scale_transfert_function[2] = GTK_RANGE(
				lookup_widget("scale_highlights"));
	}

	assert(fit->naxes[2] == 1 || fit->naxes[2] == 3);
	nb_chan = fit->naxes[2];
	ndata = fit->rx * fit->ry;

	m = gtk_range_get_value(scale_transfert_function[1]);
	lo = gtk_range_get_value(scale_transfert_function[0]);
	hi = gtk_range_get_value(scale_transfert_function[2]);
	pente = 1.0 / (hi - lo);

	undo_save_state("Processing: Histogram Transformation (mid=%.3lf, low=%.3lf, high=%.3lf)",
			m, lo, hi);

	for (chan = 0; chan < nb_chan; chan++) {
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)
		for (i = 0; i < ndata; i++) {
			double pxl = ((double) buf[chan][i] / (double) norm);
			pxl -= lo;
			pxl *= pente;
			buf[chan][i] = round_to_WORD(MTF(pxl, m) * (double) norm);
		}
	}
}

gboolean on_scale_button_release_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	set_cursor_waiting(TRUE);
	update_histo_mtf();
	set_cursor_waiting(FALSE);
	return FALSE;
}

gboolean on_scale_key_release_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	set_cursor_waiting(TRUE);
	update_histo_mtf();
	set_cursor_waiting(FALSE);
	return FALSE;
}

void on_button_histo_apply_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	apply_mtf_to_fits(&gfit);
	_init_clipped_pixels();
	update_gfit_histogram_if_needed();
	reset_curors_and_values();
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

double MTF(double x, double m) {
	double out;

	if (m == 0)
		out = 0.0;
	else if (m == 0.5)
		out = x;
	else if (m == 1)
		out = 1;
	else {
		out = ((m - 1.0) * x) / (((2.0 * m - 1.0) * x) - m);
	}
	return out;
}

void apply_mtf_to_histo(gsl_histogram *histo, double norm, double m, double lo,
		double hi) {
	gsl_histogram *mtf_histo;
	unsigned short i;
	double pente = 1.0 / (hi - lo);

	mtf_histo = gsl_histogram_alloc((size_t) norm + 1);
	gsl_histogram_set_ranges_uniform(mtf_histo, 0, norm);

#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)
	for (i = 0; i < round_to_WORD(norm); i++) {
		WORD mtf;
		double binval = gsl_histogram_get(histo, i);
		double pxl = ((double) i / norm);
		uint64_t clip[2] = { 0, 0 };

		if (i < round_to_WORD(lo * norm)) {
			pxl = lo;
			clip[0] += binval;
		} else if (i > round_to_WORD(hi * norm)) {
			pxl = hi;
			clip[1] += binval;
		}
		pxl -= lo;
		pxl *= pente;
		mtf = round_to_WORD(MTF(pxl, m) * norm);
		gsl_histogram_accumulate(mtf_histo, mtf, binval);
#pragma omp critical
		{
			clipped[0] += clip[0];
			clipped[1] += clip[1];
		}
	}
	gsl_histogram_memcpy(histo, mtf_histo);
	gsl_histogram_free(mtf_histo);
}

void update_histo_mtf() {
	double lo, hi, mid;
	unsigned int i, data = 0;
	static GtkRange *scale_transfert_function[3] = { NULL, NULL, NULL };
	static GtkAdjustment *adjustment[3] = { NULL, NULL, NULL };
	static GtkWidget *drawarea = NULL, *selarea = NULL;
	double norm = (double) gsl_histogram_bins(com.layers_hist[0]) - 1;

	if (drawarea == NULL) {
		drawarea = lookup_widget("drawingarea_histograms");
		selarea = lookup_widget("drawingarea_histograms_selection");
	}
	if (scale_transfert_function[0] == NULL) {
		scale_transfert_function[0] = GTK_RANGE(lookup_widget("scale_shadows"));
		scale_transfert_function[1] = GTK_RANGE(
				lookup_widget("scale_midtones"));
		scale_transfert_function[2] = GTK_RANGE(
				lookup_widget("scale_highlights"));
		adjustment[0] = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adj_shadows"));
		adjustment[1] = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adj_midtones"));
		adjustment[2] = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adj_highlights"));
	}

	lo = gtk_range_get_value(scale_transfert_function[0]);
	mid = gtk_range_get_value(scale_transfert_function[1]);
	hi = gtk_range_get_value(scale_transfert_function[2]);

	for (i = 0; i < 3; i++)
		gtk_adjustment_set_step_increment(adjustment[i], 0.00001);

	gtk_range_set_range(scale_transfert_function[0], 0.0, hi);
	gtk_range_set_range(scale_transfert_function[2], lo, 1.0);

	update_gfit_histogram_if_needed();	// take time. Need to be optimized

	_init_clipped_pixels();
	for (i = 0; i < gfit.naxes[2]; i++) {
		apply_mtf_to_histo(com.layers_hist[i], norm, mid, lo, hi);
	}
	data = gfit.rx * gfit.ry * gfit.naxes[2];
	_update_clipped_pixels(data);

	/* redraw the histogram */
	gtk_widget_queue_draw(drawarea);
	gtk_widget_queue_draw(selarea);
}
