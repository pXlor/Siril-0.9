#ifndef STACKING_H_
#define STACKING_H_

/* the dynamic image selection, based on various possible criteria */
typedef int (*stack_filter)(sequence *seq, int nb_img, double param);

/* the stacking method */
struct stacking_args;
typedef int (*stack_method)(struct stacking_args *args);

typedef struct normalization_coeff norm_coeff;


/* TYPE OF SIGMA CLIPPING */
typedef enum {
	NO_REJEC,
	PERCENTILE,
	SIGMA,
	SIGMEDIAN,
	WINSORIZED,
	LINEARFIT,
} rejection;

/* TYPE OF NORMALIZATION */
typedef enum {
	NO_NORM,
	ADDITIVE,
	MULTIPLICATIVE,
	ADDITIVE_SCALING,
	MULTIPLICATIVE_SCALING,
} normalization;

struct normalization_coeff {
	double *offset;
	double *mul;
	double *scale;
};

struct stacking_args {
	stack_method method;
	sequence *seq;
	stack_filter filtering_criterion;
	double filtering_parameter;
	int nb_images_to_stack; // calculated from the above, for display purposes
	int *image_indices;	// conversion between selected image indices and sequence image indices
	char description[100];
	const char *output_filename;
	gboolean output_overwrite;
	struct timeval t_start;
	int retval;
	int max_number_of_rows;	/* number of rows that can be processed simultaneously,
				   function of max memory, image size and nb_images_to_stack */
	double sig[2];		/* low and high sigma rejection */
	rejection type_of_rejection;		/* Type of rejection */
	normalization normalize;		/* Normalization */
	gboolean force_norm;		/* TRUE = force normalization */
};

void initialize_stacking_methods();

int stack_summing(struct stacking_args *args);
int stack_median(struct stacking_args *args);
int stack_mean_with_rejection(struct stacking_args *args);
int stack_addmax(struct stacking_args *args);
int stack_addmin(struct stacking_args *args);

void start_stacking();
void update_stack_interface();

#endif
