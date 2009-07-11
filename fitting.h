#ifndef INC_FITTING_H
#define INC_FITTING_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include "commath.h"

typedef struct {
	double (*f)(const float x, const float y, const gsl_vector *vars,
		const float *fpars);
	float *poss;
	int xlen;
	int *yb;
	float *fpars; // a1 and a2
} Errparams;

double fixedparab(const float x, const float y, const gsl_vector *vars,
		const float *as);
void minerror(const Errparams *pars, const size_t n, gsl_vector *vars);

#endif /* INC_FITTING_H */
