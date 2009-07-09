#ifndef INC_COMMATH_H
#define INC_COMMATH_H

#include <stdio.h> // Required for size_t data type
#include <math.h>
#include <gsl/gsl_cblas.h>

/* STRUCTURES
 * **********
 */
typedef struct {
    size_t degree;
    float *coeffs;
} Polynom;

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
//float fmax(float a, float b);
//float fmin(float a, float b);


#endif /* INC_COMMATH_H */
