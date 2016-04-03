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

#define MAIN
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <gtk/gtk.h>
#include <unistd.h>
#include <signal.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#if (defined(__APPLE__) && defined(__MACH__))
#include <stdlib.h>
#include <libproc.h>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "io/conversion.h"
#include "gui/callbacks.h"
#include "registration/registration.h"
#include "stacking/stacking.h"
#include "core/undo.h"
#include "io/single_image.h"

/* the global variables of the whole project */
cominfo com;	// the main data struct
fits gfit;	// currently loaded image
fits wfit[5];	// used for temp files, can probably be replaced by local variables
GtkBuilder *builder;	// get widget references anywhere
void initialize_scrollbars();


char *siril_sources[] = {
	"",
#if (defined(__APPLE__) && defined(__MACH__))
	"",
#endif
	PACKAGE_DATA_DIR"/",
	"/usr/share/siril/",
	"/usr/local/share/siril/"
};

void usage(const char *command) {
    printf("\nUsage:  %s [OPTIONS] [IMAGE_FILE_TO_OPEN]\n\n", command);
    puts("-i                      With init file name in argument. Start Siril.");
    puts("-f (or --format)        Print all supported image file formats (depending on the libraries you've installed)");
    puts("-v (or --version)       Print program name and version and exit");
    puts("-h (or --help)          This text");
}

void signal_handled(int s) {
	//printf("Caught signal %d\n", s);
	undo_flush();
	exit(EXIT_FAILURE);

}

int main(int argc, char *argv[]) {
	int i;
	extern char *optarg;
	extern int opterr;
	GString *homeStr;
	gchar *homepath, *csspath, *siril_path;
	char *cwd_orig = NULL;
	struct sigaction sigIntHandler;
#if (defined(__APPLE__) && defined(__MACH__))
	int ret;
	pid_t pid;
	char path[PROC_PIDPATHINFO_MAXSIZE];
#endif
	
	setenv("LC_NUMERIC", "C", 1);		// avoid possible bugs using french separator ","
	setenv("LANG", "C", 1);				// on french (or other) system, it avoids the mix between languages
	opterr = 0;
	memset(&com, 0, sizeof(struct cominf));	// needed?
	com.initfile = NULL;

	/* Caught signals */
	sigIntHandler.sa_handler = signal_handled;
	sigemptyset(&sigIntHandler.sa_mask);
	sigIntHandler.sa_flags = 0;
	sigaction(SIGINT, &sigIntHandler, NULL);

	while(1) {
		signed char c = getopt(argc, argv, "i:hfv");
		if (c=='?') {
			for(i=1; i<argc; i++) {
				if(argv[i][1] == '-') {
					if(!strcmp(argv[i], "--version"))
						c='v';
					else if(!strcmp(argv[i], "--help"))
						c='h';				
					else if(!strcmp(argv[i], "--format"))
						c='f';					
					else {
						usage(argv[0]);
						exit(EXIT_FAILURE);
					}
				}
			}
		}

		if (c==-1) break;
		switch(c) {
			case 'i':
				com.initfile = strdup(optarg);
				break;
			case 'v':
				fprintf(stdout, "%s %s\n", PACKAGE, VERSION);
				exit(EXIT_SUCCESS);
				break;
			case 'f':
				list_format_available();
				exit(EXIT_SUCCESS);
				break;
			default:
				fprintf(stderr, "unknown command line parameter `%c'\n", c);
				/* no break */
			case 'h':
				usage(argv[0]);
				exit(EXIT_SUCCESS);
		}
	}

	gtk_init (&argc, &argv);

#if (defined(__APPLE__) && defined(__MACH__))
	pid = getpid();
	ret = proc_pidpath ( pid, path, sizeof(path));
	if ( ret ) {
		int l;
		char* sp;
		sp = strrchr ( path, '/' );
		sp -= 5;
		*sp = 0;
		l = strlen ( path );
		if (( l + 10 ) < PROC_PIDPATHINFO_MAXSIZE ) {
			( void ) strcat ( path, "Resources/" );
			siril_sources[1] = path;
		}
	}
#endif

	/* try to load the glade file, from the sources defined above */
	builder = gtk_builder_new();
	i = 0;
	do {
		GError *err = NULL;
		GString *pathStr = g_string_new (siril_sources[i]);
		g_string_append(pathStr, "siril3.glade");
		gchar *path = g_string_free (pathStr, FALSE);

		if (gtk_builder_add_from_file (builder, path, &err)) {
			fprintf(stdout, "successfully loaded %s\n", siril_sources[i]);
			g_free(path);
			break;
		}
		g_free(path);
		g_error_free(err);
		i++;
	} while (i < sizeof(siril_sources)/sizeof(char *));
	if (i == sizeof(siril_sources) / sizeof(char *)) {
		fprintf(stderr, "siril3.glade was not found or contains errors, cannot render GUI. Exiting.\n");
		exit(EXIT_FAILURE);
	}
	siril_path = siril_sources[i];

	gtk_builder_connect_signals (builder, NULL);

	/* Create tags associated with the buffer for the output text. */
	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(gtk_builder_get_object(builder, "output")));
	/* Tag with weight bold and tag name "bold" . */
	gtk_text_buffer_create_tag (tbuf, "bold", "weight", PANGO_WEIGHT_BOLD, NULL);
	/* Tag with style normal */
	gtk_text_buffer_create_tag (tbuf, "normal", "weight", PANGO_WEIGHT_NORMAL, NULL);
	/* Couleur Tags */
	gtk_text_buffer_create_tag (tbuf, "red", "foreground", "#e72828", NULL);
	gtk_text_buffer_create_tag (tbuf, "salmon", "foreground", "#ff9898", NULL);
	gtk_text_buffer_create_tag (tbuf, "green", "foreground", "#01b301", NULL);
	gtk_text_buffer_create_tag (tbuf, "blue", "foreground", "#7a7af8", NULL);
	gtk_text_buffer_create_tag (tbuf, "plum", "foreground", "#8e4585", NULL);
	
	siril_log_color_message("Welcome to %s v%s\n", "bold", PACKAGE, VERSION);

	/* initialize converters (utilities used for different image types importing) */
	initialize_converters();

	/* initializing internal structures with widgets (drawing areas) */
	com.vport[RED_VPORT] = lookup_widget("drawingarear");
	com.vport[GREEN_VPORT] = lookup_widget("drawingareag");
	com.vport[BLUE_VPORT] = lookup_widget("drawingareab");
	com.vport[RGB_VPORT] = lookup_widget("drawingareargb");
	com.preview_area[0] = lookup_widget("drawingarea_preview1");
	com.preview_area[1] = lookup_widget("drawingarea_preview2");
	com.cvport = RED_VPORT;
	com.show_excluded = TRUE;
	com.selected_star = -1;
	com.star_is_seqdata = FALSE;
	com.stars = NULL;
	com.uniq = NULL;
	com.grad = NULL;
	com.grad_boxes_drawn = TRUE;
	com.color = NORMAL_COLOR;
	for (i=0; i<MAXVPORT; i++)
		com.buf_is_dirty[i] = TRUE;
	initialize_remap();
	memset(&com.selection, 0, sizeof(rectangle));
	memset(com.layers_hist, 0, sizeof(com.layers_hist));
	initialize_scrollbars();
	init_mouse();

	/* Keybord Shortcuts */
	initialize_shortcuts();

	/* Fill the "About" window */
	fill_about_dialog();

	/* Select combo boxes that trigger some text display or other things */
	gtk_combo_box_set_active(GTK_COMBO_BOX(gtk_builder_get_object(builder, "comboboxstack_methods")), 0);
	gtk_combo_box_set_active(GTK_COMBO_BOX(gtk_builder_get_object(builder, "comboboxstacksel")), 0);

	/* initialize the com struct and zoom level */
	com.drawn = FALSE;
	com.sliders = MINMAX;
	com.zoom_value = ZOOM_DEFAULT;
	zoomcombo_update_display_for_zoom();
	
	/* initialize comboboxs of extraction background */
	GtkComboBox *order = GTK_COMBO_BOX(gtk_builder_get_object(builder, "combo_polyorder"));
	gtk_combo_box_set_active(order, POLY_4);
	GtkComboBox *grad_inter = GTK_COMBO_BOX(gtk_builder_get_object(builder, "combobox_gradient_inter"));
	gtk_combo_box_set_active(grad_inter, 0);

	/* initialize sequence-related stuff */
	initialize_sequence(&com.seq, TRUE);
	adjust_sellabel();

	/* set default CWD and load init file */
	com.wd = malloc(PATH_MAX + 1);// PATH_MAX may not be available on all systems
	if (!getcwd(com.wd, PATH_MAX)) {
		free(com.wd);
		com.wd = NULL;
	} else
		cwd_orig = strdup(com.wd);

	if (checkinitfile()) {
		siril_log_message("could not load or create settings file in ~/.siril, exiting.\n");
		exit(1);
	}
	
	/* load the css sheet for general style */
	homepath = getenv("HOME");
	homeStr = g_string_new (homepath);
	g_string_append(homeStr, "/.siril");
	csspath = g_string_free (homeStr, FALSE);

	load_css_style_sheet (siril_path, csspath);

	g_free(csspath);

	/* initialize registration methods */
	initialize_registration_methods();

	/* initialize stacking methods */
	initialize_stacking_methods();

	/* initialization of the binning parameters */
	GtkComboBox *binning = GTK_COMBO_BOX(gtk_builder_get_object(builder, "combobinning"));
	gtk_combo_box_set_active(binning, 0);

	/* initialization of swap dir chooser */
	GtkFileChooser *swap_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_swap"));
	GtkLabel *label = GTK_LABEL(lookup_widget("label_swap_dir"));

	if (com.swap_dir && com.swap_dir[0] != '\0')
		gtk_file_chooser_set_filename (swap_dir, com.swap_dir);
	else
		gtk_file_chooser_set_filename (swap_dir, g_get_tmp_dir());
	gtk_label_set_text (label, com.swap_dir);

	/* initialization of default FITS extension */
	GtkComboBox *box = GTK_COMBO_BOX(lookup_widget("combobox_ext"));
	gtk_combo_box_set_active_id(box, com.ext);

#ifdef HAVE_LIBRAW
	set_GUI_LIBRAW();
#endif
	
	/* Get CPU number and set the number of threads */
	siril_log_message("Parallel processing %s: Using %d logical processor%s",
#ifdef _OPENMP
			"enabled", com.max_thread = omp_get_num_procs(), "s.\n"
#else
			"disabled", com.max_thread = 1, ".\n"
#endif
			);
	update_spinCPU(com.max_thread);

	/* start Siril */
	update_used_memory();

	if (argv[optind] != NULL) {
		if (cwd_orig)
			changedir(cwd_orig);
		open_single_image(argv[optind]);
		char *newpath = extract_path(argv[optind]);
		changedir(newpath);
		free(newpath);
	}
	if (cwd_orig)
		free(cwd_orig);
	gtk_main();

	/* quit Siril */
	undo_flush();
	return 0;
}

void initialize_scrollbars() {
	int i;
	char *vport_names[] = { "r", "g", "b", "rgb" };
	char window_name[32];

	for (i=0; i<sizeof(vport_names)/sizeof(char *); i++) {
		sprintf(window_name, "scrolledwindow%s", vport_names[i]);
		GtkScrolledWindow *win = GTK_SCROLLED_WINDOW(gtk_builder_get_object(builder, window_name));
		com.hadj[i] = gtk_scrolled_window_get_hadjustment(win);
		g_signal_connect(com.hadj[i], "value-changed",
				G_CALLBACK(scrollbars_hadjustment_changed_handler), NULL);
		com.vadj[i] = gtk_scrolled_window_get_vadjustment(win);
		g_signal_connect(com.vadj[i], "value-changed",
				G_CALLBACK(scrollbars_vadjustment_changed_handler), NULL);
	}
}
