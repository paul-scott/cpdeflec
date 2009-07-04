#ifndef INC_COMMATH_H
#define INC_COMMATH_H

#include <stdio.h> // Required for size_t data type
#include <gsl/gsl_cblas.h>

/* STRUCTURES
 * **********
 */
typedef struct {
    size_t degree;
    float coeffs[];
} Polynom;

/* FUNCTIONS
 * *********
 */
float getvalpoly(Polynom *pol, float val);
float fdot(const float *x, const float *y);
float fnorm(const float *x);
float fscale(const float alpha, float *x);


#endif /* INC_COMMATH_H */
