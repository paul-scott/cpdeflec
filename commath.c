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

#include "commath.h"

const double PI = 3.14159265358979;


double polyget(Polynom *pol, double val)
{
	double out = 0.0;
	for (int i=0; i<(pol->degree+1); i++) {
		out = out + *(pol->coeffs+i)*pow(val, (double) i);
	}
	return out;
}

// Could potentially use macros to implement these:
double dot(const double *x, const double *y)
{
	return cblas_ddot(3, x, 1, y, 1);
}

double norm(const double *x)
{
	return cblas_dnrm2(3, x, 1);
}

void scale(const double alpha, double *x)
{
	return cblas_dscal(3, alpha, x, 1);
}

void axpy(const double alpha, const double *x, double *y)
{
	cblas_daxpy(3, alpha, x, 1, y, 1);
}

void cross(const double *a, const double *b, double *c)
{
	*(c) = (*(a+1))*(*(b+2)) - (*(a+2))*(*(b+1));
	*(c+1) = (*(a+2))*(*(b)) - (*(a))*(*(b+2));
	*(c+2) = (*(a))*(*(b+1)) - (*(a+1))*(*(b));
}

void matxvec(const double *M, const double *x, double *y)
{
	// Note that these functions could first copy data so that data could
	// be copied back to first variable.
	*(y) = dot(M,x);
	*(y+1) = dot(M+1*3,x);
	*(y+2) = dot(M+2*3,x);
}

void matxmat(const double *m1, const double *m2, double *m3)
{
	*(m3) = (*(m1))*(*(m2)) + (*(m1+1))*(*(m2+3)) + (*(m1+2))*(*(m2+6));
	*(m3+1) = (*(m1))*(*(m2+1)) + (*(m1+1))*(*(m2+4)) + (*(m1+2))*(*(m2+7));
	*(m3+2) = (*(m1))*(*(m2+2)) + (*(m1+1))*(*(m2+5)) + (*(m1+2))*(*(m2+8));
	*(m3+3) = (*(m1+3))*(*(m2)) + (*(m1+4))*(*(m2+3)) + (*(m1+5))*(*(m2+6));
	*(m3+4) = (*(m1+3))*(*(m2+1)) + (*(m1+4))*(*(m2+4)) + (*(m1+5))*(*(m2+7));
	*(m3+5) = (*(m1+3))*(*(m2+2)) + (*(m1+4))*(*(m2+5)) + (*(m1+5))*(*(m2+8));
	*(m3+6) = (*(m1+6))*(*(m2)) + (*(m1+7))*(*(m2+3)) + (*(m1+8))*(*(m2+6));
	*(m3+7) = (*(m1+6))*(*(m2+1)) + (*(m1+7))*(*(m2+4)) + (*(m1+8))*(*(m2+7));
	*(m3+8) = (*(m1+6))*(*(m2+2)) + (*(m1+7))*(*(m2+5)) + (*(m1+8))*(*(m2+8));
}

int imax(int a, int b)
{
	return a > b ? a : b;
}

int imin(int a, int b)
{
	return a < b ? a : b;
}

void rotmatxyz(const double xa, const double ya, const double za, double *rm)
{
	// Creates rotation vector from applying rx then ry then rz.
	double *m1 = malloc(9*sizeof(*m1));
	double *m2 = malloc(9*sizeof(*m2));

	rotmatx(xa, rm);
	rotmaty(ya, m1);
	matxmat(m1, rm, m2);
	rotmatz(za, m1);
	matxmat(m1, m2, rm);

	free(m1);
	m1 = NULL;
	free(m2);
	m2 = NULL;
}
void rotmatzyx(const double za, const double ya, const double xa, double *rm)
{
	// Creates rotation vector from applying rz then ry then rx.
	double *m1 = malloc(9*sizeof(*m1));
	double *m2 = malloc(9*sizeof(*m2));

	rotmatz(za, rm);
	rotmaty(ya, m1);
	matxmat(m1, rm, m2);
	rotmatx(xa, m1);
	matxmat(m1, m2, rm);

	free(m1);
	m1 = NULL;
	free(m2);
	m2 = NULL;
}

void rotmatx(const double ang, double *rm)
{
	*(rm) = 1.0;
	*(rm+1) = 0.0;
	*(rm+2) = 0.0;
	*(rm+3) = 0.0;
	*(rm+4) = cos(ang);
	*(rm+5) = sin(ang);
	*(rm+6) = 0.0;
	*(rm+7) = -sin(ang);
	*(rm+8) = cos(ang);
}

void rotmaty(const double ang, double *rm)
{
	*(rm) = cos(ang);
	*(rm+1) = 0.0;
	*(rm+2) = -sin(ang);
	*(rm+3) = 0.0;
	*(rm+4) = 1.0;
	*(rm+5) = 0.0;
	*(rm+6) = sin(ang);
	*(rm+7) = 0.0;
	*(rm+8) = cos(ang);
}

void rotmatz(const double ang, double *rm)
{
	*(rm) = cos(ang);
	*(rm+1) = sin(ang);
	*(rm+2) = 0.0;
	*(rm+3) = -sin(ang);
	*(rm+4) = cos(ang);
	*(rm+5) = 0.0;
	*(rm+6) = 0.0;
	*(rm+7) = 0.0;
	*(rm+8) = 1.0;
}
