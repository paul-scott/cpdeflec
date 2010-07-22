#ifndef INC_FITTING_H
#define INC_FITTING_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include "commath.h"

typedef struct {
	double *poss;
	int xlen;
	int *yb;
	double *fpars; // a1 and a2
} Errparams;

double sphesqerr(const gsl_vector *vars, void *params);
double parabsqerr(const gsl_vector *vars, void *params);
void minerror(double (*f)(const gsl_vector *va, void *params),
		const Errparams *pars, const size_t n, gsl_vector *vars,
		const gsl_vector *stepsize);
double paraboloid(const double x, const double y, const double a,
		const double b);
double sphere(const double x, const double y, const double rad);
void sphereslope(const double x, const double y, const double rad, double *nrm);
void parabslope(const double x, const double y, const double f1,
		const double f2, double *nrm);

#endif /* INC_FITTING_H */
