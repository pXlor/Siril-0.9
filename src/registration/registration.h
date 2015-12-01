#ifndef _REGISTRATION_H_
#define _REGISTRATION_H_

#include "core/siril.h"

struct registration_args;
typedef int (*registration_function)(struct registration_args *);

/* arguments passed to registration functions */
struct registration_args {
	registration_function func;	// the registration function
	sequence *seq;			// the sequence to register
	gboolean process_all_frames;	// all frames of the sequence (opposite of selected frames)
	rectangle selection;		// the selection rectangle
	int layer;			// layer of images on which the registration is computed
	struct timeval t_start;		// start time of func
	int retval;			// retval of func
	gboolean run_in_thread;		// true if the registration was run in a thread
	GtkEntry *entry;
	const gchar *text;
};


typedef enum { REQUIRES_NO_SELECTION,	// selection is not used
	REQUIRES_ANY_SELECTION,		// selection can be of any size and shape
	REQUIRES_SQUARED_SELECTION	// selection needs to be square-shaped
} selection_type;


/* used to register a registration method */
struct registration_method {
	const char *name;
	registration_function method_ptr;
	selection_type sel;
};


struct registration_method *new_reg_method(char *name, registration_function f, selection_type s); // for compositing
void	initialize_registration_methods();
struct registration_method * get_selected_registration_method();
int	register_shift_dft(struct registration_args *args);
int	register_shift_fwhm(struct registration_args *args);
int register_star_alignment(struct registration_args *args);
int	get_registration_layer(sequence *seq);
void	update_reg_interface(gboolean dont_change_reg_radio);
void	get_the_registration_area(struct registration_args *reg_args, struct registration_method *method); // for compositing


/* mouse behaviour */

typedef enum {
	MOUSE_ACTION_NONE,
	MOUSE_ACTION_SELECT_REG_AREA,
	MOUSE_ACTION_SELECT_PREVIEW1,
	MOUSE_ACTION_SELECT_PREVIEW2,
} mouse_status_enum;

extern mouse_status_enum mouse_status;	// defined in registration_preview.c

#endif