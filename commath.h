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

#ifndef INC_COMMATH_H
#define INC_COMMATH_H

#include <stdio.h> // Required for size_t data type
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_cblas.h>

/* STRUCTURES
 * **********
 */
typedef struct {
    size_t degree;
    double *coeffs;
} Polynom;

extern const double PI;

/* FUNCTIONS
 * *********
 */
double polyget(Polynom *pol, double val);
double dot(const double *x, const double *y);
double norm(const double *x);
void scale(const double alpha, double *x);
void axpy(const double alpha, const double *x, double *y);
void cross(const double *a, const double *b, double *c);
void matxvec(const double *M, const double *x, double *y);
int imax(int a, int b);
int imin(int a, int b);
void matxmat(const double *m1, const double *m2, double *m3);
void rotmatxyz(const double xa, const double ya, const double za, double *rm);
void rotmatzyx(const double za, const double ya, const double xa, double *rm);
void rotmatx(const double ang, double *rm);
void rotmaty(const double ang, double *rm);
void rotmatz(const double ang, double *rm);


#endif /* INC_COMMATH_H */
