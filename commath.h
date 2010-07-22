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
