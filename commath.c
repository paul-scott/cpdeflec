#include "commath.h"

const double PI = 3.14159265358979;


float polyget(Polynom *pol, float val) {
	float out = 0.0f;
	for (int i=0; i<(pol->degree+1); i++) {
		out = out + *(pol->coeffs+i)*powf(val, (float) i);
	}
	return out;
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

void fmatxmat(const float *m1, const float *m2, float *m3) {
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

int imax(int a, int b) {
	return a > b ? a : b;
}

int imin(int a, int b) {
	return a < b ? a : b;
}

void rotmatxyz(const float xa, const float ya, const float za, float *rm) {
	// Creates rotation vector from applying rx then ry then rz.
	float *m1 = malloc(9*sizeof(*m1));
	float *m2 = malloc(9*sizeof(*m2));

	rotmatx(xa, rm);
	rotmaty(ya, m1);
	fmatxmat(m1, rm, m2);
	rotmatz(za, m1);
	fmatxmat(m1, m2, rm);

	free(m1);
	m1 = NULL;
	free(m2);
	m2 = NULL;
}

void rotmatx(const float ang, float *rm) {
	*(rm) = 1.0f;
	*(rm+1) = 0.0f;
	*(rm+2) = 0.0f;
	*(rm+3) = 0.0f;
	*(rm+4) = cosf(ang);
	*(rm+5) = sinf(ang);
	*(rm+6) = 0.0f;
	*(rm+7) = -sinf(ang);
	*(rm+8) = cosf(ang);
}

void rotmaty(const float ang, float *rm) {
	*(rm) = cosf(ang);
	*(rm+1) = 0.0f;
	*(rm+2) = -sinf(ang);
	*(rm+3) = 0.0f;
	*(rm+4) = 1.0f;
	*(rm+5) = 0.0f;
	*(rm+6) = sinf(ang);
	*(rm+7) = 0.0f;
	*(rm+8) = cosf(ang);
}

void rotmatz(const float ang, float *rm) {
	*(rm) = cosf(ang);
	*(rm+1) = sinf(ang);
	*(rm+2) = 0.0f;
	*(rm+3) = -sinf(ang);
	*(rm+4) = cosf(ang);
	*(rm+5) = 0.0f;
	*(rm+6) = 0.0f;
	*(rm+7) = 0.0f;
	*(rm+8) = 1.0f;
}
