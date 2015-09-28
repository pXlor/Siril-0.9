#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "siril.h"

typedef struct
{
	size_t row, col;
	size_t box;
	size_t NbBoxes;
	size_t box_per_row;
	size_t box_per_col;
	double tolerance;
	double deviation;
	double unbalance;
} Im2sub;

gsl_matrix *Bkg_1color(gsl_vector *, Im2sub *);
int extract_background(fits *, fits*, int);
double poly_1(gsl_vector *, int, int);
double poly_2(gsl_vector *, int, int);
double poly_3(gsl_vector *, int, int);
double poly_4(gsl_vector *, int, int);
void extract_BKG(fits *);

void fill_boxes_list();
void add_box_to_list(gradient grad, int);
