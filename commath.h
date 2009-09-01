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
    float *coeffs;
} Polynom;

extern const double PI;

/* FUNCTIONS
 * *********
 */
float polyget(Polynom *pol, float val);
float fdot(const float *x, const float *y);
float fnorm(const float *x);
void fscale(const float alpha, float *x);
void faxpy(const float alpha, const float *x, float *y);
void fcross(const float *a, const float *b, float *c);
void fmatxvec(const float *M, const float *x, float *y);
int imax(int a, int b);
int imin(int a, int b);
void fmatxmat(const float *m1, const float *m2, float *m3);
void rotmatxyz(const float xa, const float ya, const float za, float *rm);
void rotmatx(const float ang, float *rm);
void rotmaty(const float ang, float *rm);
void rotmatz(const float ang, float *rm);


#endif /* INC_COMMATH_H */
