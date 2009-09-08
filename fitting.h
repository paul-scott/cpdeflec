#ifndef INC_FITTING_H
#define INC_FITTING_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include "commath.h"

typedef struct {
	float *poss;
	int xlen;
	int *yb;
	float *fpars; // a1 and a2
} Errparams;

double sphesqerr(const gsl_vector *vars, void *params);
double parabsqerr(const gsl_vector *vars, void *params);
void minerror(double (*f)(const gsl_vector *va, void *params),
		const Errparams *pars, const size_t n, gsl_vector *vars,
		const gsl_vector *stepsize);
double paraboloid(const double x, const double y, const double a,
		const double b);
double sphere(const double x, const double y, const double rad);
void sphereslope(const float x, const float y, const float rad, float *nrm);
void parabslope(const float x, const float y, const float f1, const float f2,
		float *nrm);

#endif /* INC_FITTING_H */
