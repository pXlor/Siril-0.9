#include <assert.h>
#include <string.h>
#include <glib.h>

#include "siril.h"
#include "processing.h"
#include "proto.h"
#include "gui/callbacks.h"
#include "io/ser.h"

// called in start_in_new_thread only
// works in parallel if the arg->parallel is TRUE for FITS or SER sequences
gpointer generic_sequence_worker(gpointer p) {
	struct generic_seq_args *args = (struct generic_seq_args *) p;
	struct timeval t_end;
	float nb_framesf, progress = 0.f; // 0 to nb_framesf, for progress
	int frame;	// the current frame, sequence index
	int current; // number of processed frames so far
	int abort = 0; // variable for breaking out of loop
	GString *desc; // temporary string description for logs
	gchar *msg; // final string description for logs

	assert(args);
	assert(args->seq);
	assert(args->image_hook);
	set_progress_bar_data(NULL, PROGRESS_RESET);
	if (args->nb_filtered_images > 0)
		nb_framesf = (float)args->nb_filtered_images;
	else 	nb_framesf = (float)args->seq->number;
	args->retval = 0;

	if (args->prepare_hook && args->prepare_hook(args)) {
		siril_log_message("Preparing sequence processing failed.\n");
		args->retval = 1;
		goto the_end;
	}

	/* Output print of algorithm description */
	desc = g_string_new(args->description);
	if (desc) {
		desc = g_string_append(desc, ": processing...\n");
		msg = g_string_free(desc, FALSE);
		siril_log_color_message(msg, "red");
		g_free(msg);
	}

	current = 0;

#pragma omp parallel for num_threads(com.max_thread) private(frame) schedule(static) \
	if(args->parallel && ((args->seq->type == SEQ_REGULAR && fits_is_reentrant()) || args->seq->type == SEQ_SER))
	for (frame = 0; frame < args->seq->number; frame++) {
		if (!abort) {
			fits fit;
			char filename[256], msg[256];

			if (!get_thread_run()) {
				abort = 1;
				continue;
			}
			if (args->filtering_criterion
					&& !args->filtering_criterion(args->seq, frame,
							args->filtering_parameter))
				continue;

			if (!seq_get_image_filename(args->seq, frame, filename)) {
				abort = 1;
				continue;
			}

			snprintf(msg, 256, "%s. Processing image %d (%s)", args->description, frame, filename);
			progress = (float) (args->nb_filtered_images <= 0 ? frame : current);
			set_progress_bar_data(msg, progress / nb_framesf);

			memset(&fit, 0, sizeof(fits));
			if (seq_read_frame(args->seq, frame, &fit)) {
				abort = 1;
				clearfits(&fit);
				continue;
			}

			if (args->image_hook(args, frame, current, &fit)) {
				abort = 1;
				clearfits(&fit);
				continue;
			}

			int retval;
			if (args->save_hook)
				retval = args->save_hook(args, frame, current, &fit);
			else retval = generic_save(args, frame, current, &fit);
			if (retval) {
				abort = 1;
				clearfits(&fit);
				continue;
			}

			clearfits(&fit);

#pragma omp atomic
			current++;
		}
	}

	if (abort) {
		set_progress_bar_data("Sequence processing failed. Check the log.", PROGRESS_RESET);
		siril_log_message("Sequence processing failed.\n");
		args->retval = abort;
	}
	else {
		set_progress_bar_data("Sequence processing succeeded.", PROGRESS_RESET);
		siril_log_message("Sequence processing succeeded.\n");
		gettimeofday(&t_end, NULL);
		show_time(args->t_start, t_end);
	}

the_end:
	if (args->finalize_hook && args->finalize_hook(args)) {
		siril_log_message("Finalizing sequence processing failed.\n");
		args->retval = 1;
	}

	if (args->idle_function)
		gdk_threads_add_idle(args->idle_function, args);
	else gdk_threads_add_idle(end_generic_sequence, args);
	return NULL;
}

// defaut idle function (in GTK main thread) to run at the end of the generic sequence processing
gboolean end_generic_sequence(gpointer p) {
	struct generic_seq_args *args = (struct generic_seq_args *) p;

	if (args->load_new_sequence && args->new_seq_prefix && !args->retval) {
		char *seqname = malloc(strlen(args->new_seq_prefix) + strlen(args->seq->seqname) + 5);
		sprintf(seqname, "%s%s.seq", args->new_seq_prefix, args->seq->seqname);
		check_seq(0);
		update_sequences_list(seqname);
		free(seqname);
	}
	
	return end_generic(p);
}

int ser_prepare_hook(struct generic_seq_args *args) {
	char dest[256];
	const char *ptr;

	if (args->force_ser_output || args->seq->type == SEQ_SER) {
		ptr = strrchr(args->seq->seqname, '/');
		if (ptr)
			snprintf(dest, 255, "%s%s.ser", args->new_seq_prefix, ptr + 1);
		else snprintf(dest, 255, "%s%s.ser", args->new_seq_prefix, args->seq->seqname);

		args->new_ser = malloc(sizeof(struct ser_struct));
		return ser_create_file(dest, args->new_ser, TRUE, args->seq->ser_file);
	}

	return 0;
}

int ser_finalize_hook(struct generic_seq_args *args) {
	int retval = 0;
	if (args->force_ser_output || args->seq->type == SEQ_SER) {
		retval = ser_write_and_close(args->new_ser);
		free(args->new_ser);
	}
	return retval;
}

int generic_save(struct generic_seq_args *args, int input_index, int processed_index, fits *fit) {
	char dest[256];
	if (args->force_ser_output || args->seq->type == SEQ_SER) {
		snprintf(dest, 256, "%s%s.ser", args->new_seq_prefix, args->seq->seqname);
		return ser_write_frame_from_fit(args->new_ser, fit, input_index);
	} else {
		snprintf(dest, 256, "%s%s%05d%s", args->new_seq_prefix,
				args->seq->seqname, processed_index, com.ext);
		return savefits(dest, fit);
	}
}

/*****************************************************************************
 *      P R O C E S S I N G      T H R E A D      M A N A G E M E N T        *
 ****************************************************************************/

// This function is reentrant
void start_in_new_thread(gpointer (*f)(gpointer p), gpointer p) {
	g_mutex_lock(&com.mutex);
	if (com.run_thread || com.thread != NULL) {
		fprintf(stderr, "The processing thread is busy, stop it first.\n");
		g_mutex_unlock(&com.mutex);
		return;
	}

	com.run_thread = TRUE;
	g_mutex_unlock(&com.mutex);
	com.thread = g_thread_new("processing", f, p);
}

void stop_processing_thread() {
	if (com.thread == NULL) {
		fprintf(stderr,
				"The processing thread is not running, cannot stop it.\n");
		return;
	}

	set_thread_run(FALSE);

	g_thread_join(com.thread);
	com.thread = NULL;
}

void set_thread_run(gboolean b) {
	g_mutex_lock(&com.mutex);
	com.run_thread = b;
	g_mutex_unlock(&com.mutex);
}

gboolean get_thread_run() {
	gboolean retval;
	g_mutex_lock(&com.mutex);
	retval = com.run_thread;
	g_mutex_unlock(&com.mutex);
	return retval;
}

/* should be called in a threaded function if nothing special has to be done at the end.
 * gdk_threads_add_idle(end_generic, NULL);
 */
gboolean end_generic(gpointer arg) {
	stop_processing_thread();
	update_used_memory();
	set_cursor_waiting(FALSE);
	return FALSE;
}

void on_processes_button_cancel_clicked(GtkButton *button, gpointer user_data) {
	if (com.thread != NULL)
		siril_log_color_message("Process aborted by user\n", "red");
	stop_processing_thread();
}

int seq_filter_all(sequence *seq, int nb_img, double any) {
	return 1;
}

int seq_filter_included(sequence *seq, int nb_img, double any) {
	return (seq->imgparam[nb_img].incl);
}


