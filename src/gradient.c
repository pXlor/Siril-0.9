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
#include <assert.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <string.h>

#include "callbacks.h"
#include "gradient.h"
#include "siril.h"
#include "proto.h"


#define NPARAM_POLY4 15		// Number of parameters used with 4rd order
#define NPARAM_POLY3 10		// Number of parameters used with 3rd order
#define NPARAM_POLY2 6		// Number of parameters used with 2nd order
#define NPARAM_POLY1 3		// Number of parameters used with 1nd order

gsl_matrix *Bkg_1color(gsl_vector *MatR, Im2sub *Info){
	size_t i, j, k, n;
	size_t row, col, inc_row, inc_col;
	size_t inc=0;
	double chisq, pixel_value, tmp_row, tmp_col;
	double tolerance = Info->tolerance;
	double deviation = Info->deviation;
	double unbalance = Info->unbalance;
	gsl_matrix *J, *cov;
	gsl_vector *y, *w, *c;
	size_t box = Info->box, midbox = box * 0.5;
	size_t box_per_row = Info->box_per_row;
	size_t box_per_col = Info->box_per_col;
	size_t height = Info->row;
	size_t width = Info->col;	
	
	assert(box > 0);
	
	n = box_per_row * box_per_col;
	if (com.grad) free(com.grad);
	com.grad = malloc(n * sizeof (gradient));
	
	GtkComboBox *combo_order = GTK_COMBO_BOX(gtk_builder_get_object(builder, "combo_polyorder"));
	poly_order order = gtk_combo_box_get_active(combo_order);
	
	if (box_per_row < (int)order + 1 || box_per_col < (int)order + 1) {
		return NULL;
	}
	
	gsl_vector *Vec_row = gsl_vector_alloc(box_per_col);
	gsl_vector *Vec_col = gsl_vector_alloc(box_per_row);
	
	tmp_row=(double)midbox - 1.0;
	tmp_col=(double)midbox - 1.0;

	for (i=0; i<box_per_col; i++){
		gsl_vector_set(Vec_row, i, tmp_row);
		tmp_row += (double)((height - 2 * midbox) / (box_per_col - 1));
	}
	for (i=0; i<box_per_row; i++){
		gsl_vector_set(Vec_col, i, tmp_col);
		tmp_col += (double)((width - 2 * midbox) / (box_per_row - 1));
	}

	gsl_vector *mailles_row = gsl_vector_alloc(n); //line coordinates of the box center
	gsl_vector *mailles_col = gsl_vector_alloc(n); //column coordinates of the box center
	gsl_vector *maillesR = gsl_vector_alloc(n); //mean of the pixel intensity in the box. Pixels above threshold are rejected of the mean calcul.
	
	for (row=0; row<box_per_col; row++) {
		size_t start_row = round(gsl_vector_get(Vec_row, row) - midbox + 1);
		for (col=0; col<box_per_row; col++){
			double *data_box = calloc(box*box, sizeof(double));
			size_t start_col = round(gsl_vector_get(Vec_col, col) - midbox + 1);
			for (inc_row=0, k=0; inc_row<box; inc_row++) {
				for (inc_col=0; inc_col<box; inc_col++) {
					data_box[k++] = gsl_vector_get(MatR, ((size_t)start_row+inc_row) * width + start_col + inc_col);
				}
			}
			double sigma = gsl_stats_sd(data_box, 1, k);
			gsl_sort(data_box, 1, k);
			double median = gsl_stats_median_from_sorted_data(data_box, 1, k);
			
			for (inc_row=0, k=0; inc_row<box; inc_row++) {
				for (inc_col=0; inc_col<box; inc_col++) {
					double current_pixel = gsl_vector_get(MatR, (size_t)(start_row + inc_row) * width + start_col + inc_col);
					if (current_pixel > (tolerance * sigma + median)) {
						gsl_vector_set(MatR, (size_t)(start_row + inc_row) * width + start_col + inc_col, median);
					}
					data_box[k++] = gsl_vector_get(MatR,(size_t)(start_row + inc_row) * width + start_col + inc_col);
				}
			}
			gsl_sort(data_box, 1, k);
			double value = gsl_stats_median_from_sorted_data(data_box, 1, k);
			
			gsl_vector_set(maillesR, inc, value);
								
			com.grad[inc].centre.x = gsl_vector_get(Vec_col, col) + midbox;
			com.grad[inc].centre.y = height - gsl_vector_get(Vec_row, row) + midbox;	
			gsl_vector_set(mailles_row, inc, gsl_vector_get(Vec_row, row));
			gsl_vector_set(mailles_col, inc, gsl_vector_get(Vec_col, col));
			inc++;
			free(data_box);
		}
	}
	
	double *data = calloc(n, sizeof(double));
	for (i=0; i<n; i++)
		data[i] = gsl_vector_get(maillesR, i);
	
	gsl_sort(data, 1, n);
	double median = gsl_stats_median_from_sorted_data(data, 1, n);
	double sigma = gsl_stats_sd(data, 1, n);
	
	for (i=0; i<n; i++) {
		double pixel = gsl_vector_get(maillesR, i);
		if (((pixel - median)/sigma > deviation) || ((median - pixel)/sigma > (deviation * unbalance)))
			gsl_vector_set(maillesR, i, -1);
		com.grad[i].boxvalue = gsl_vector_get(maillesR, i);
	}
	free(data);	

	
	int VarParam;
	switch (order) {
		case POLY_1:
			VarParam = NPARAM_POLY1;
			break;
		case POLY_2:
			VarParam = NPARAM_POLY2;
			break;
		case POLY_3:
			VarParam = NPARAM_POLY3;
			break;
		case POLY_4:
		default:
			VarParam = NPARAM_POLY4;
	}

	// J is the Jacobian
	// y contains data (pixel intensity)
	J = gsl_matrix_alloc (n, VarParam);
	y = gsl_vector_alloc (n);
	w = gsl_vector_alloc (n);
	c = gsl_vector_alloc (VarParam);
	cov = gsl_matrix_alloc (VarParam, VarParam);
	
	for (inc = 0; inc < n; inc++) {
		tmp_col = gsl_vector_get(mailles_col, inc);
		tmp_row = gsl_vector_get(mailles_row, inc);
		pixel_value  = gsl_vector_get(maillesR, inc);
		//here, it is a bit sketchy in the sense that if there is not value to report in a box (because the threshold is too
		//low for example, then I just skip the initialization of J and y. gsl automatically discard the non assigned values
		//during the minimization. I tested it with Matlab and it works fine. The results agree.
		if (pixel_value < 0) continue;
		
		gsl_matrix_set (J, inc, 0, 1.0);
		gsl_matrix_set (J, inc, 1, tmp_col);
		gsl_matrix_set (J, inc, 2, tmp_row);
		
		if (order != POLY_1) {
			gsl_matrix_set (J, inc, 3, tmp_col * tmp_col);
			gsl_matrix_set (J, inc, 4, tmp_col * tmp_row);
			gsl_matrix_set (J, inc, 5, tmp_row * tmp_row);
		}
		
		if (order == POLY_3 || order == POLY_4) {
			gsl_matrix_set (J, inc, 6, tmp_col * tmp_col * tmp_col);
			gsl_matrix_set (J, inc, 7, tmp_col * tmp_col * tmp_row);
			gsl_matrix_set (J, inc, 8, tmp_col * tmp_row * tmp_row);
			gsl_matrix_set (J, inc, 9, tmp_row * tmp_row * tmp_row);
		}
		
		if (order == POLY_4) {
			gsl_matrix_set (J, inc, 10, tmp_col * tmp_col * tmp_col * tmp_col);
			gsl_matrix_set (J, inc, 11, tmp_col * tmp_col * tmp_col * tmp_row);
			gsl_matrix_set (J, inc, 12, tmp_col * tmp_col * tmp_row * tmp_row);
			gsl_matrix_set (J, inc, 13, tmp_col * tmp_row * tmp_row * tmp_row);
			gsl_matrix_set (J, inc, 14, tmp_row * tmp_row * tmp_row * tmp_row);
		}
		
		gsl_vector_set (y, inc, pixel_value);
		gsl_vector_set (w, inc, 1.0);
	}

	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, VarParam);
	gsl_multifit_wlinear (J, w, y, c, cov, &chisq, work);
	gsl_multifit_linear_free (work);
	gsl_matrix_free	(J);
	gsl_vector_free	(y);
	gsl_vector_free	(w);
	gsl_vector_free	(maillesR);
	gsl_vector_free	(mailles_row);
	gsl_vector_free	(mailles_col);

//Calcul of the background with the same dimension that the input matrix.
	gsl_matrix *Bkg = gsl_matrix_alloc(height,width);

	for (i=0; i<height; i++) {
		for (j=0; j<width; j++) {
			switch (order) {
				case POLY_1:
					pixel_value = poly_1(c, j, i);
					break;
				case POLY_2:
					pixel_value = poly_2(c, j, i);
					break;
				case POLY_3:
					pixel_value = poly_3(c, j, i);
					break;
				case POLY_4:
				default:		//default should not happend
					pixel_value = poly_4(c, j, i);
			}
			gsl_matrix_set(Bkg, i, j, pixel_value);
		}
	}
	return Bkg;
}

int build_double_array(double *dbl_array, const gsl_vector *vec_src){
	size_t incdata;
	for (incdata=0; incdata < vec_src->size; incdata++)	{
		dbl_array[incdata] = gsl_vector_get(vec_src,incdata);
	}
	return 0;
}

int find_indices(gsl_vector *vec_indices, gsl_vector *vec_src, double elem){
	size_t inc;
	size_t count = 0;
	for (inc=0; inc < vec_src->size; inc++)	{
		//printf("vec = %lg, elem = %lg\n",gsl_vector_get(vec_src,inc),elem);
		if (gsl_vector_get(vec_src, inc) == elem)		{
			gsl_vector_set(vec_indices, count, (double) inc);
			count++;
		}
	}
	return (int) count;
}

int cpy_all_but_one(gsl_vector *vec_dest, const gsl_vector *vec_src, const size_t pos){
	size_t inc;
	size_t count=0;
	double tmp;
	for (inc=0; inc < vec_src->size; inc++)	{
		if (inc == pos) continue;
		tmp = gsl_vector_get(vec_src, inc);
		gsl_vector_set(vec_dest, count, tmp);
		count++;
	}
	return 0;
}


int NNAveraging(gsl_vector *maillesX, int NbBoxesc, int NbBoxesl,
			int oneptx, int onepty, int posptrem, gsl_vector *xrem, gsl_vector *yrem, int neighbour) {
	gsl_vector *xremallorcut;
	gsl_vector *yremallorcut;
	gsl_vector *xremtrans;
	gsl_vector *yremtrans;
	gsl_vector *reducedxrem;
	gsl_vector *reducedyrem;
	gsl_vector *xotherpt2rem;
	gsl_vector *PtsSelectedx;
	gsl_vector *PtsSelectedy;
	gsl_vector *IntPtsSel;
	PtsSelectedx = gsl_vector_alloc(NbBoxesc * NbBoxesl);
	PtsSelectedy = gsl_vector_alloc(NbBoxesc * NbBoxesl);
	int incptrem, incx, incy;
	int increduced;
	size_t incfin;
	
	if(oneptx == 0) {
		xremallorcut = gsl_vector_alloc(xrem->size);
		yremallorcut = gsl_vector_alloc(yrem->size);
		gsl_vector_memcpy(xremallorcut, xrem);
		gsl_vector_memcpy(yremallorcut, yrem);
	}
	else {
		xremallorcut = gsl_vector_alloc(1);
		yremallorcut = gsl_vector_alloc(1);
		gsl_vector_set(xremallorcut,0,oneptx);
		gsl_vector_set(yremallorcut,0,onepty);
	}
	
	for(incptrem = 0; incptrem<(int) xremallorcut->size; incptrem++) {
		
		if(oneptx==0)
			posptrem = incptrem;
		int ptremx = gsl_vector_get(xremallorcut,(size_t) incptrem);
		int ptremy = gsl_vector_get(yremallorcut,(size_t) incptrem);
		int inc = 0;
		int startloopx, endloopx, startloopy, endloopy;
		
		startloopx = ptremx-neighbour-1;
		endloopx = ptremx+neighbour-1;
		startloopy = ptremy-neighbour-1;
		endloopy = ptremy+neighbour-1;
		
		if (startloopx < 0) startloopx = 0;
		if (endloopx >= NbBoxesc) endloopx = NbBoxesc-1;
		if (startloopy < 0) startloopy = 0;
		if (endloopy >= NbBoxesl) endloopy = NbBoxesl-1;

		for (incx=startloopx; incx<=endloopx; incx++) {
			for (incy=startloopy; incy<=endloopy; incy++) {
				if (incx == (ptremx-1) && incy == (ptremy-1)) continue;
				//~ flag=0;
				xremtrans = gsl_vector_alloc(xrem->size);
				yremtrans = gsl_vector_alloc(yrem->size);
				reducedxrem = gsl_vector_alloc(xrem->size - 1);
				reducedyrem = gsl_vector_alloc(yrem->size - 1);
				xotherpt2rem = gsl_vector_alloc(reducedxrem->size);
				gsl_vector_memcpy(xremtrans,xrem);
				gsl_vector_memcpy(yremtrans,yrem);
				gsl_vector_add_constant(xremtrans,-1);
				gsl_vector_add_constant(yremtrans,-1);
				cpy_all_but_one(reducedxrem,xremtrans,(size_t) posptrem);
				cpy_all_but_one(reducedyrem,yremtrans,(size_t) posptrem);
				int count = find_indices(xotherpt2rem, reducedxrem, (double) incx);
				if (count > 0) {
					for (increduced=0; increduced < count; increduced++) {
						size_t posxrem = (size_t) gsl_vector_get(xotherpt2rem,(size_t) increduced);
						if ((int) gsl_vector_get(reducedyrem,posxrem) == incy) {
							//~ flag = 1;
							continue;
						}
					}
					//~ if (flag==1) {
						//~ flag=0;
						//~ continue;
					//~ }
				}
				gsl_vector_set(PtsSelectedx, inc, incx);
				gsl_vector_set(PtsSelectedy, inc, incy);
				inc++;
				gsl_vector_free(xremtrans);
				gsl_vector_free(yremtrans);
				gsl_vector_free(reducedxrem);
				gsl_vector_free(reducedyrem);
				gsl_vector_free(xotherpt2rem);
			}
		}
	
		if (inc>0) {
			IntPtsSel = gsl_vector_alloc((size_t) inc);
			for (incfin=0; incfin < inc; incfin++) {
				size_t posX = (size_t) gsl_vector_get(PtsSelectedx, incfin);
				size_t posY = (size_t) gsl_vector_get(PtsSelectedy, incfin);
				gsl_vector_set(IntPtsSel,incfin,gsl_vector_get(maillesX, posX + posY * NbBoxesc));
			}
			double data[inc];
			build_double_array(data, IntPtsSel);
			gsl_vector_set(maillesX, (size_t) ptremx-1 + (ptremy-1) * NbBoxesc, gsl_stats_mean(data, 1, inc));
			gsl_vector_free(IntPtsSel);
		}
		if (inc<4) {
			NNAveraging(maillesX, NbBoxesc, NbBoxesl, ptremx, ptremy, posptrem, xrem, yrem, neighbour+1);
		}
	}
	
	gsl_vector_free(xremallorcut);
	gsl_vector_free(yremallorcut);
	gsl_vector_free(PtsSelectedx);
	gsl_vector_free(PtsSelectedy);
	 
	return 0;
}

gsl_matrix *Bkg_1color_spline(gsl_vector *MatR, Im2sub *Info) {
	int i, n, k;
	double tmp, tmpl, tmpc;
	double tolerance = Info->tolerance;
	double deviation = Info->deviation;
	double unbalance = Info->unbalance;
	int box = Info->box, midbox = box * 0.5;
	size_t box_per_row = Info->box_per_row;
	size_t box_per_col = Info->box_per_col;
	size_t height = Info->row;
	size_t width = Info->col;
	size_t row, col, inc_col, inc_row, inc;		
	gsl_vector *Vec_row = gsl_vector_alloc(box_per_col);
	gsl_vector *Vec_col = gsl_vector_alloc(box_per_row);
	
	tmpl=(double)midbox - 1;
	tmpc=(double)midbox - 1;
	
	n = box_per_row * box_per_col;
	if (com.grad) free(com.grad);
	com.grad = malloc(n * sizeof (gradient));
	
	for (i=0; i<box_per_col; i++) {
		gsl_vector_set(Vec_row, i, tmpl);
		tmpl += (double) (height - 2 * midbox) / (box_per_col - 1);
	}
	for (i=0; i<box_per_row; i++) {
		gsl_vector_set(Vec_col, i, tmpc);
		tmpc += (double) (width - 2 * midbox) / (box_per_row - 1);
	}
	
	gsl_vector *mailles_row = gsl_vector_alloc(n);
	gsl_vector *mailles_col = gsl_vector_alloc(n);
	gsl_vector *maillesR = gsl_vector_alloc(n);
	
	for (row=0, inc=0; row<box_per_col; row++) {
		size_t start_row = round(gsl_vector_get(Vec_row, row) - midbox + 1);
		for (col=0; col<box_per_row; col++){
			double *data_box = calloc(box * box, sizeof(double));
			size_t start_col = round(gsl_vector_get(Vec_col, col) - midbox + 1);
			for (inc_row=0, k=0; inc_row<box; inc_row++) {
				for (inc_col=0; inc_col<box; inc_col++) {
					data_box[k++] = gsl_vector_get(MatR, ((size_t)start_row + inc_row) * width + start_col + inc_col);
				}
			}
			double sigma = gsl_stats_sd(data_box, 1, k);
			gsl_sort(data_box, 1, k);
			double median = gsl_stats_median_from_sorted_data(data_box, 1, k);
			
			for (inc_row=0, k=0; inc_row<box; inc_row++) {
				for (inc_col=0; inc_col<box; inc_col++)	{
					double current_pixel = gsl_vector_get(MatR, (size_t)(start_row + inc_row) * width + start_col + inc_col);
					if (current_pixel < (tolerance * sigma + median))
						data_box[k++] = current_pixel;
				}
			}
			gsl_sort(data_box, 1, k);
			median = gsl_stats_median_from_sorted_data(data_box, 1, k);
			
			gsl_vector_set(maillesR, inc, median);
			gsl_vector_set(mailles_row, inc, gsl_vector_get(Vec_row, row));
			gsl_vector_set(mailles_col, inc, gsl_vector_get(Vec_col, col));
			com.grad[inc].centre.x = gsl_vector_get(Vec_col, col) + midbox;
			com.grad[inc].centre.y = height - gsl_vector_get(Vec_row, row) + midbox;	
			com.grad[inc].boxvalue = gsl_vector_get(maillesR, inc);
			inc++;
			free(data_box);
		}
	}
	
	double *data = calloc(n, sizeof(double));
	for (i=0; i<n; i++)
		data[i] = gsl_vector_get(maillesR, i);

	double sigma = gsl_stats_sd(data, 1, n);
	gsl_sort(data, 1, n);
	double median = gsl_stats_median_from_sorted_data(data, 1, n);
	
	size_t *x_tab = calloc(n, sizeof(size_t));
	size_t *y_tab = calloc(n, sizeof(size_t));
	
	size_t box_removed=0;
	for (row=0; row<box_per_col; row++) {
		for (col=0; col<box_per_row; col++){
			double pixel = gsl_vector_get(maillesR, col + (box_per_row * row));
			if (((pixel - median)/sigma > deviation) || ((median - pixel)/sigma > (deviation * unbalance))) {
				x_tab[box_removed] = col + 1;
				y_tab[box_removed] = box_per_col - row;
				com.grad[col + (box_per_row * row)].boxvalue = -1.0;
				box_removed++;
			}
			//~ printf("%lu %lu: %g, %g\n", col + 1, row + 1, pixel, (pixel - median) / sigma);
		}
	}
	free(data);	
	fill_boxes_list(); 
	
	if (box_removed > 1) {
		gsl_vector *xrem = gsl_vector_alloc(box_removed);
		gsl_vector *yrem = gsl_vector_alloc(box_removed);
		
		for (i = 0; i < box_removed; i++) {
			gsl_vector_set(xrem, i, x_tab[i]);
			gsl_vector_set(yrem, i, y_tab[i]);
			//~ printf("%lu et %lu\n", x_tab[i], y_tab[i]);
		}
		
		int oneptx = 0;
		int onepty = 0;
		int posptrem = 0;
		int neighbour = 1;
		
		NNAveraging(maillesR, (int) box_per_row, (int) box_per_col, oneptx, onepty, posptrem, xrem, yrem, neighbour);
		
		gsl_vector_free(xrem);
		gsl_vector_free(yrem);
	}
	free(x_tab);
	free(y_tab);
		
	/* smoothing maillesR */
	gsl_vector *maillesTmp = gsl_vector_alloc(n);
	gsl_vector_memcpy(maillesTmp, maillesR);
	for (row=1; row<box_per_col - 1; row++) {
		for (col=1; col<box_per_row - 1; col++){
			double tmp = gsl_vector_get(maillesR, col + (box_per_row * row));
			tmp += gsl_vector_get(maillesR, col + (box_per_row * (row - 1 )));
			tmp += gsl_vector_get(maillesR, col + 1 + (box_per_row * (row - 1)));
			tmp += gsl_vector_get(maillesR, col + 1 + (box_per_row * row));
			tmp += gsl_vector_get(maillesR, col + 1 + (box_per_row * (row + 1)));
			tmp += gsl_vector_get(maillesR, col + (box_per_row * (row + 1)));
			tmp += gsl_vector_get(maillesR, col - 1 + (box_per_row * (row + 1)));
			tmp += gsl_vector_get(maillesR, col - 1 + (box_per_row * row));
			tmp += gsl_vector_get(maillesR, col - 1 + (box_per_row * (row - 1)));
			gsl_vector_set(maillesTmp, col + (box_per_row * row), tmp / 9.0);
		}
	}
	gsl_vector_memcpy(maillesR, maillesTmp);
	gsl_vector_free(maillesTmp);
						
	gsl_matrix *Bkg = gsl_matrix_alloc(height, width);
	double f00, f01, f10, f11, tmpcBeg, tmpcEnd, tmplBeg, tmplEnd, ValMeshX, ValMeshY;
	size_t sizeMeshX, sizeMeshY, incMeshX, incMeshY;

	for (row = 0; row < box_per_col-1; row++) {
		for (col = 0; col < box_per_row-1; col++) {
			tmplBeg = round(gsl_vector_get(mailles_col, (row + 1) * box_per_row + col));
			tmplEnd = round(gsl_vector_get(mailles_col, (row + 1) * box_per_row + col + 1));
			tmpcBeg = round(gsl_vector_get(mailles_row, row * box_per_row + col));
			tmpcEnd = round(gsl_vector_get(mailles_row, (row + 1) * box_per_row + col + 1));
								
			f00  = gsl_vector_get(maillesR, row * box_per_row + col);
			f10  = gsl_vector_get(maillesR, row * box_per_row + col + 1);
			f01  = gsl_vector_get(maillesR, (row + 1) * box_per_row + col);
			f11  = gsl_vector_get(maillesR, (row + 1) * box_per_row + col + 1);
			
			sizeMeshY = tmplEnd - tmplBeg + 1;
			sizeMeshX = tmpcEnd - tmpcBeg + 1;
			
			for (incMeshY = 0; incMeshY < sizeMeshY; incMeshY++){
				ValMeshY = incMeshY / (tmplEnd - tmplBeg);
				for (incMeshX = 0; incMeshX < sizeMeshX; incMeshX++){
					ValMeshX = incMeshX / (tmpcEnd - tmpcBeg);
					tmp = (1 - ValMeshX) * (1 - ValMeshY) * f00;
					tmp+= ValMeshX * (1 - ValMeshY) * f01;
					tmp+= (1 - ValMeshX) * ValMeshY * f10;
					tmp+= ValMeshX * ValMeshY * f11;
					gsl_matrix_set(Bkg, incMeshX + tmpcBeg, incMeshY + tmplBeg, tmp);
				}
			}
        }
    }

    /* Compute Edges */

	double begEdges = 0;
	double begEdgesY = width - midbox;
	double endEdges = round(gsl_vector_get(Vec_row, 1)) + midbox;

	for (row=0; row<box_per_col - 1; row++) {
		tmpcEnd = endEdges;
		tmplEnd = midbox - 1;
	
		f00  = gsl_vector_get(maillesR, row * box_per_row);
		f10  = gsl_vector_get(maillesR, row * box_per_row);
		f01  = gsl_vector_get(maillesR, (row + 1) * box_per_row);
		f11  = gsl_vector_get(maillesR, (row + 1) * box_per_row);
		
		sizeMeshY = tmplEnd + 1;
		sizeMeshX = tmpcEnd + 1;
		for (incMeshX = begEdges; incMeshX < sizeMeshX; incMeshX++)	{
			ValMeshX = (incMeshX - begEdges) / (endEdges - begEdges);
			for (incMeshY = 0; incMeshY < sizeMeshY; incMeshY++) {
				ValMeshY = (double)incMeshY / (sizeMeshY - 1);
				tmp = (1 - ValMeshX) * (1 - ValMeshY) * f00;
				tmp+= ValMeshX * (1 - ValMeshY) * f01;
				tmp+= (1 - ValMeshX) * ValMeshY * f10;
				tmp+= ValMeshX * ValMeshY * f11;
				gsl_matrix_set(Bkg , incMeshX, incMeshY, tmp);
			}
		}

		f00  = gsl_vector_get(maillesR, (row + 1) * (box_per_row) - 1);
		f10  = gsl_vector_get(maillesR, (row + 1) * (box_per_row) - 1);
		f01  = gsl_vector_get(maillesR, (row + 1) * (box_per_row) + box_per_row - 1);
		f11  = gsl_vector_get(maillesR, (row + 1) * (box_per_row) + box_per_row - 1);
		
		for (incMeshX = begEdges; incMeshX < sizeMeshX; incMeshX++)	{
			ValMeshX = (incMeshX-begEdges) / (endEdges - begEdges);
			for (incMeshY = begEdgesY; incMeshY < width; incMeshY++){
				ValMeshY = (double)(incMeshY-begEdgesY) / (midbox - 1);
				tmp = (1 - ValMeshX) * (1 - ValMeshY) * f00;
				tmp+= ValMeshX * (1 - ValMeshY) * f01;
				tmp+= (1 - ValMeshX) * ValMeshY * f10;
				tmp+= ValMeshX * ValMeshY * f11;
				gsl_matrix_set(Bkg , incMeshX, incMeshY, tmp);
			}
		}
	
		if (row < (box_per_col - 2)){
			begEdges = endEdges;
			endEdges = round(gsl_vector_get(Vec_row, row + 2)) + midbox;
		}
	}
	
	begEdges = 0;
	double begEdgesX = height - midbox;

	endEdges = round(gsl_vector_get(Vec_col, 1)) + midbox;
	for (col=0; col<box_per_row - 1; col++) {
		tmplEnd = endEdges;
		tmpcEnd = midbox - 1;
	
		f00  = gsl_vector_get(maillesR, col);
		f10  = gsl_vector_get(maillesR, col + 1);
		f01  = gsl_vector_get(maillesR, col);
		f11  = gsl_vector_get(maillesR, col + 1);
		sizeMeshY = tmplEnd + 1;
		sizeMeshX = tmpcEnd + 1;
		for (incMeshY = begEdges; incMeshY < sizeMeshY; incMeshY++)	{
			ValMeshY = (incMeshY-begEdges) / (endEdges - begEdges);
			for (incMeshX = 0; incMeshX < sizeMeshX; incMeshX++){
				ValMeshX = (double) incMeshX / (sizeMeshX - 1);
				tmp = (1 - ValMeshX) * (1 - ValMeshY) * f00;
				tmp+= ValMeshX * (1 - ValMeshY) * f01;
				tmp+= (1 - ValMeshX) * ValMeshY * f10;
				tmp+= ValMeshX * ValMeshY * f11;
				gsl_matrix_set(Bkg , incMeshX, incMeshY, tmp);
			}
		}
	
		f00  = gsl_vector_get(maillesR, (col + box_per_row * (box_per_col - 1)));
		f10  = gsl_vector_get(maillesR, (col + 1 + box_per_row * (box_per_col - 1)));
		f01  = gsl_vector_get(maillesR, (col + box_per_row * (box_per_col - 1)));
		f11  = gsl_vector_get(maillesR, (col + 1 + box_per_row * (box_per_col - 1)));
	
		for (incMeshY = begEdges; incMeshY < sizeMeshY; incMeshY++){
			ValMeshY = (incMeshY-begEdges) / (endEdges - begEdges);
			for (incMeshX = begEdgesX; incMeshX < height; incMeshX++){
				ValMeshX = (double) (incMeshX - begEdgesX) / (midbox - 1);
				tmp = (1 - ValMeshX) * (1 - ValMeshY) * f00;
				tmp+= ValMeshX * (1 - ValMeshY) * f01;
				tmp+= (1 - ValMeshX) * ValMeshY * f10;
				tmp+= ValMeshX * ValMeshY * f11;
				gsl_matrix_set(Bkg , incMeshX, incMeshY, tmp);
			}
		}
	
		if (col < (box_per_row - 2)) {
			begEdges = endEdges;
			endEdges = round(gsl_vector_get(Vec_col, col + 2)) + midbox;
		}
	}
	return Bkg;
}


//C contains background function
#define C(i) (gsl_vector_get(c,(i)))

double poly_4(gsl_vector * c, int x, int y) {
	double value = C(0)*1 + C(1)*x + C(2)*y + C(3)*x*x + C(4)*y*x + C(5)*y*y +
							C(6)*x*x*x + C(7)*x*x*y + C(8)*x*y*y + C(9)*y*y*y +
							C(10)*x*x*x*x + C(11)*x*x*x*y + C(12)*x*x*y*y + C(13)*x*y*y*y + C(14)*y*y*y*y;
	if (value<0.) value=0.;
	return value;
}

double poly_3(gsl_vector * c, int x, int y) {
	double value = C(0)*1 + C(1)*x + C(2)*y + C(3)*x*x + C(4)*y*x + C(5)*y*y +
							C(6)*x*x*x + C(7)*x*x*y + C(8)*x*y*y + C(9)*y*y*y;
	if (value<0.) value=0.;
	return value;
}

double poly_2(gsl_vector * c, int x, int y) {
	double value = C(0)*1 + C(1)*x + C(2)*y + C(3)*x*x + C(4)*y*x + C(5)*y*y;
	if (value<0.) value=0.;	
	return value;
}

double poly_1(gsl_vector * c, int x, int y) {
	double value = C(0)*1 + C(1)*x + C(2)*y;
	if (value<0.) value=0.;	
	return value;
}

int extract_background(fits *imgfit, fits *bkgfit, int layer) {
	Im2sub *Info = malloc(sizeof(Im2sub));
	size_t ndata, i, j;
	WORD *buf = imgfit->pdata[layer];
	
	ndata=(size_t)(imgfit->rx*imgfit->ry);
	gsl_vector *MatR = gsl_vector_alloc(ndata);
	gsl_matrix *Bkg;
	
	Info->row		= (size_t)imgfit->ry; 
	Info->col		= (size_t)imgfit->rx; 
	Info->box		= (size_t)gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_bkg_sizebox"))) * 2;
	int interval	= (size_t)gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_bkg_Box_sep")));     
	double tolerance = (double)gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_bkg_tolerance")));
	double deviation = (double)gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_bkg_deviation")));
	double unbalance = (double)gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_bkg_unbalance")));
	int inter_type = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combobox_gradient_inter")));
	
	Info->box_per_row	= (size_t)((double)Info->col / ((double)Info->box + interval - 1));
	Info->box_per_col	= (size_t)((double)Info->row / ((double)Info->box + interval - 1));

	Info->tolerance = tolerance;
	Info->deviation = deviation;
	Info->unbalance = unbalance;
		
	for (i=0; i<ndata; i++) {
		gsl_vector_set(MatR,i,(double)buf[i]);
	}
		
	if (inter_type==1) {
		Info->box_per_row	= (size_t)gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_bkg_Box_s_p_row")));
		Info->box_per_col	= (size_t)(((double)Info->box_per_row / imgfit->rx) * imgfit->ry);
		
		/* Fill Global Variables */
		com.grad_nb_boxes = Info->box_per_col * Info->box_per_row;
		com.grad_size_boxes = Info->box;
			
		if (!(Bkg = Bkg_1color_spline(MatR, Info))) {
			return 1;
		} 
	}
	else {
		/* Fill Global Variables */
		com.grad_nb_boxes = Info->box_per_col * Info->box_per_row;
		com.grad_size_boxes = Info->box;
		
		if (!(Bkg = Bkg_1color(MatR, Info))) {
			return 1;							// not enough samples
		} 
	}
	gsl_vector_free(MatR);
	
	if (imgfit->naxes[2]>1)
		copyfits(imgfit, bkgfit, CP_ALLOC | CP_FORMAT | CP_EXPAND, layer);
	else
		copyfits(imgfit, bkgfit, CP_ALLOC | CP_FORMAT | CP_COPYA, 0);
	
	WORD *tbuf = bkgfit->pdata[layer];
	for (i=0; i<Info->row; i++) {
		for (j=0; j<Info->col; j++)
			*tbuf++=(WORD)gsl_matrix_get(Bkg, i, j);
	}
		
	siril_log_message("Channel #%d: background extraction done.\n", layer);
	gsl_matrix_free(Bkg);
	return 0;
}

void extract_BKG(fits *fit){
	int i;
	struct timeval t_start, t_end;

	siril_log_color_message("Background extraction: processing...\n", "red");	
	gettimeofday (&t_start, NULL);	
		
	for (i=0; i<com.uniq->nb_layers; i++)
		if (extract_background(&gfit, fit, i)) {
			siril_log_message("Insufficient background samples.\n");
			return;
		}
	gtk_widget_set_sensitive(lookup_widget("frame_bkg_tools"), TRUE);
	gtk_widget_set_sensitive(lookup_widget("button_bkg_correct"), TRUE);
	
	gettimeofday (&t_end, NULL);	
	show_time(t_start, t_end);
	
}

/******* Boxes list ******************/
static GtkListStore *liststore_bkg = NULL;

enum {
	COLUMN_INDEX,			// gint
	COLUMN_SELECTED,		// gboolean
	N_COLUMNS
};


void get_boxes_list_store() {
	if (liststore_bkg == NULL)
		liststore_bkg = GTK_LIST_STORE(gtk_builder_get_object(builder, "liststore_bkg"));
}

void fill_boxes_list() {
	int i=0;
	get_boxes_list_store();
	
	gtk_list_store_clear(liststore_bkg);
	while (i < com.grad_nb_boxes) {
		add_box_to_list(com.grad[i], i);
		i++;
	}
}
 
void add_box_to_list(gradient grad, int index) {
	static GtkTreeSelection *boxes_tree = NULL;
	GtkTreeIter iter;
	gboolean is_selected = grad.boxvalue == -1.0 ? FALSE : TRUE;

	if (!boxes_tree)
		boxes_tree = GTK_TREE_SELECTION(gtk_builder_get_object(builder, "treeview-selection4"));
			
	gtk_list_store_append (liststore_bkg, &iter);
	gtk_list_store_set (liststore_bkg, &iter,
			COLUMN_INDEX, index + 1,
			COLUMN_SELECTED, is_selected,
			-1);
}

void update_selection(gchar *path, gboolean new_value) {
	GtkTreeIter iter;
	get_boxes_list_store();
	gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(liststore_bkg), &iter, path);
	gtk_list_store_set(liststore_bkg, &iter, COLUMN_SELECTED, new_value, -1);
}

void on_cellrenderertext11_toggled(GtkCellRendererToggle *cell_renderer,
		gchar *path, gpointer user_data) {
	gint *index = gtk_tree_path_get_indices(gtk_tree_path_new_from_string(path));
	if (!index) return;
	if (index[0] >= com.grad_nb_boxes) return;
	fprintf(stdout, "toggle selection index = %d\n", index[0]);

	//~ update_selection(path, !com.seq.imgparam[index[0]].incl);

	//~ com.seq.imgparam[index[0]].incl = !com.seq.imgparam[index[0]].incl;
	//~ if (com.seq.imgparam[index[0]].incl)
		//~ com.seq.selnum++;
	//~ else 	com.seq.selnum--;
	redraw(com.cvport, REMAP_NONE);
}

