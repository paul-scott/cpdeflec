#include "commath.h"


float getvalpoly(Polynom *pol, float val) {

}

// Could potentially use macros to implement these:
float fdot(const float *x, const float *y) {
	return cblas_sdot(3, x, 1, y, 1);
}

float fnorm(const float *x) {
	return cblas_snrm2(3, x, 1);
}

void fscale(const float alpha, float *x) {
	return cblas_sscal(3, alpha, x, 1);
}

void faxpy(const float alpha, const float *x, float *y) {
	cblas_saxpy(3, alpha, x, 1, y, 1);
}

void fcross(const float *a, const float *b, float *c) {
	*(c) = (*(a+1))*(*(b+2)) - (*(a+2))*(*(b+1));
	*(c+1) = (*(a+2))*(*(b)) - (*(a))*(*(b+2));
	*(c+2) = (*(a))*(*(b+1)) - (*(a+1))*(*(b));
}

void fmatxvec(const float *M, const float *x, float *y) {
	*(y) = fdot(M,x);
	*(y+1) = fdot(M+1*3,x);
	*(y+2) = fdot(M+2*3,x);
}
