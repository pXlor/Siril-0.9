#ifndef GRAD_H
#define GRAD_H

typedef struct {
	size_t row, col;
	size_t box;
	size_t NbBoxes;
	size_t box_per_row;
	size_t box_per_col;
	double tolerance;
	double deviation;
	double unbalance;
} Im2sub;

void grad_background_extraction(fits *);

#endif
