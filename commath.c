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

float fscale(const float alpha, float *x) {
	return cblas_sscal(3, alpha, x, 1);
}
