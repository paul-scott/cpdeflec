/*
 *  Copyright (C) Paul Scott 2011
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, see <http://www.gnu.org/licenses/>.
 */

#ifndef INC_FITTING_H
#define INC_FITTING_H

#include <stddef.h>
#include <gsl/gsl_vector.h>

typedef struct {
	double *poss;
	int xlen;
	int *yb;
	double *fpars; // a1 and a2
} Errparams;

double sphesqerr(const gsl_vector *vars, void *params);
double parabsqerr(const gsl_vector *vars, void *params);
void minerror(double (*f)(const gsl_vector *va, void *params),
		const Errparams *pars, size_t n, gsl_vector *vars,
		const gsl_vector *stepsize);
double paraboloid(double x, double y, double a, double b);
double sphere(double x, double y, double rad);
void sphereslope(double x, double y, double rad, double *nrm);
void parabslope(double x, double y, double f1, double f2, double *nrm);

#endif /* INC_FITTING_H */
