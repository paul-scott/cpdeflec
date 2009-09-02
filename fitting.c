#include "fitting.h"

double sphesqerr(const gsl_vector *vars, void *params) {
	// vars = xoff, yoff, zoff, rad;
	const Errparams *par = (Errparams *) params;
	float *poss = par->poss;
	int xlen = par->xlen;
	int *yb = par->yb;
	double err = 0.0f;

	for (int x=0; x<xlen; x=x+xlen/50) {
		for (int y=*(yb+x*2); y<=*(yb+x*2+1); y=y+xlen/50) {
			err = err + pow((double) (*(poss+(y*xlen+x)*3+2) -
						gsl_vector_get(vars, 2) + gsl_vector_get(vars, 3) -
						sphere(*(poss+(y*xlen+x)*3) - gsl_vector_get(vars, 0),
							*(poss+(y*xlen+x)*3+1) - gsl_vector_get(vars, 1),
							gsl_vector_get(vars, 3))),2.0);
		}
	}
	return err;
}

double parabsqerr(const gsl_vector *vars, void *params) {
	// vars = xoff, yoff, zoff, xrot, yrot, zrot, xa, ya;
	const Errparams *par = (Errparams *) params;
	float *poss = par->poss;
	int xlen = par->xlen;
	int *yb = par->yb;
	double err = 0.0f;
	float *rotm = malloc(9*sizeof(*rotm));
	float tvec1[3];
	float tvec2[3];
	// Work out rotations matrix.
	rotmatxyz((float) gsl_vector_get(vars, 3), (float) gsl_vector_get(vars, 4),
			(float) gsl_vector_get(vars, 5), rotm);

	for (int x=0; x<xlen; x=x+xlen/50) {
		for (int y=*(yb+x*2); y<=*(yb+x*2+1); y=y+xlen/50) {
			// Wait need to subtract offset from vec...
			tvec1[0] = (float) (*(poss+(y*xlen+x)*3) - gsl_vector_get(vars, 0));
			tvec1[1] = (float) (*(poss+(y*xlen+x)*3+1) -
				gsl_vector_get(vars, 1));
			tvec1[2] = (float) (*(poss+(y*xlen+x)*3+2) -
				gsl_vector_get(vars, 2));
			// Apply rotations.
			fmatxvec(rotm, tvec1, tvec2);
			// + sign since paraboloid needs to be upside down.
			err = err + pow((double) (tvec2[2] +
					paraboloid((double) tvec2[0], (double) tvec2[1],
					gsl_vector_get(vars, 6),
					gsl_vector_get(vars, 7))),2.0);
		}
	}
	free(rotm);
	rotm = NULL;
	return err;
}


double paraboloid(const double x, const double y, const double a,
		const double b) {
	// a, b is 2*foclen.
	return x*x/(2*a) + y*y/(2*b);
}

double sphere(const double x, const double y, const double rad) {
	// Need to make sure that region dosen't drift outside radius of circle.
	return sqrt(rad*rad - x*x -	y*y);
}

void sphereslope(const float x, const float y, const float rad, float *nrm) {
	// Normal is equal to sphere centre minus sphere location.
	*nrm = -x;
	*(nrm+1) = -y;
	*(nrm+2) = -sqrtf(rad*rad - x*x - y*y);
	fscale(1.0f/fnorm(nrm), nrm);
}

void minerror(double (*f)(const gsl_vector *va, void *params),
		const Errparams *pars, const size_t n, gsl_vector *vars,
		const gsl_vector *stepsize) {
	gsl_multimin_fminimizer *minzer;
	gsl_multimin_function multifunc;
	// Allocate memory for minimizer and set to nelder-mead simplex.
	minzer = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex,
			n);
	if (minzer == NULL) {
		printf("Cannot allocate memory for curve fitting...");
		exit(1);
	}

	// Set up multifunc.
	multifunc.n = n;
	multifunc.f = f;
	multifunc.params = (void *) pars;


	gsl_multimin_fminimizer_set(minzer, &multifunc, vars, stepsize);

	int status;
	size_t iter = 0;
	double size;

	do {
		status = gsl_multimin_fminimizer_iterate(minzer);

		if (status) {
			printf("Iterating failed...\n");
			break;
		}

		size = gsl_multimin_fminimizer_size(minzer);
		status = gsl_multimin_test_size(size, 1e-2);

		//printf("%5d: ", iter);
		for (int i=0; i<n; i++) {
			printf("%10.3e ", gsl_vector_get(minzer->x, i));
		}
		//printf("f() = %7.3f size = %.3f\n", minzer->fval, size);
		iter++;
	} while (status == GSL_CONTINUE && iter < 2000);

	for (int i=0; i<n; i++) {
		gsl_vector_set(vars, i, gsl_vector_get(minzer->x, i));
	}

	gsl_multimin_fminimizer_free(minzer);
	minzer = NULL;
}
