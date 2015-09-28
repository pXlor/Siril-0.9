#ifndef FWHM_H_
#define FWHM_H_
#include "siril.h"
#include "star_finder.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>

struct fwhm_struct {
	double B; /* average sky background value */
	double A; /* amplitude */
	double x0, y0; /* coordinates of the peak */
	double sx, sy; /* Size of the fitted function on the x and y axis in PSF coordinates */
	double fwhmx, fwhmy; /* FWHM in x and y axis */
	double angle; /* angle of the axis x,y with respect to the image's */
	double mag; /* magnitude of the star : this parameter is not fitted but calculated with the vector G and the parameter B */
	double xpos, ypos; /* position of the star in the image, not set by Minimization */
	double rmse; /* RMSE of the minimization */

	/* uncertainties */
	double B_err;
	double A_err;
	double x_err, y_err;
	double sx_err, sy_err;
	double ang_err;
	int layer;
	char* units;
};

struct PSF_data {
	size_t n;
	double *y;
	double *sigma;
	size_t NbRows;
	size_t NbCols;
	double rmse;
};

double get_fwhm(fits *, int, double *);
fitted_PSF *get_Minimisation(fits *, int, rectangle *);
fitted_PSF *global_Minimisation(gsl_matrix *, double, int, gboolean);
gsl_vector* MaxInz(gsl_matrix*, double);
int Gaussian_f(const gsl_vector *, void *, gsl_vector *);
int Gaussian_df(const gsl_vector *, void *, gsl_matrix *);
int Gaussian_fdf(const gsl_vector *, void *, gsl_vector *, gsl_matrix *);
int Gaussian_f_An(const gsl_vector *, void *, gsl_vector *);
int Gaussian_df_An(const gsl_vector *, void *, gsl_matrix *);
int Gaussian_fdf_An(const gsl_vector *, void *, gsl_vector *, gsl_matrix *);
fitted_PSF *minimiz_no_angle(gsl_matrix *, double, int);
fitted_PSF *minimiz_angle(gsl_matrix *, fitted_PSF *);
void DisplayResult(fitted_PSF *, rectangle *);
double Get_Magnitude(gsl_matrix *, double);
void update_units(fits*, fitted_PSF**);

#endif
