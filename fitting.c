#include "fitting.h"

double sphesqerr(const gsl_vector *vars, void *params) {
	// vars = xshft, yshft, zshft, rad;
	const Errparams *par = (Errparams *) params;
	float *poss = par->poss;
	int xlen = par->xlen;
	int *yb = par->yb;
	double err = 0.0;
	double tvec[3];

	for (int x=0; x<xlen; x=x+xlen/50) {
		for (int y=*(yb+x*2); y<=*(yb+x*2+1); y=y+xlen/50) {
			// Wait need to add shift to vec...
			tvec[0] = (double) (*(poss+(y*xlen+x)*3) +
					gsl_vector_get(vars, 0));
			tvec[1] = (double) (*(poss+(y*xlen+x)*3+1) +
					gsl_vector_get(vars, 1));
			tvec[2] = (double) (*(poss+(y*xlen+x)*3+2) +
					gsl_vector_get(vars, 2));

			err = err + pow((double) (tvec[2] + gsl_vector_get(vars, 3) -
						sphere(tvec[0], tvec[1], gsl_vector_get(vars, 3))),2.0);
		}
	}
	return err;
}

double parabsqerr(const gsl_vector *vars, void *params) {
	// vars = xshft, yshft, zshft, xrot, yrot, zrot, f1, f2;
	const Errparams *par = (Errparams *) params;
	float *poss = par->poss;
	int xlen = par->xlen;
	int *yb = par->yb;
	double err = 0.0;
	float *rotm = malloc(9*sizeof(*rotm));
	float tvec1[3];
	float tvec2[3];
	// Work out rotations matrix.
	rotmatxyz((float) gsl_vector_get(vars, 3), (float) gsl_vector_get(vars, 4),
			(float) gsl_vector_get(vars, 5), rotm);

	for (int x=0; x<xlen; x=x+xlen/50) {
		for (int y=*(yb+x*2); y<=*(yb+x*2+1); y=y+xlen/50) {
			// Wait need to add shift to vec...
			tvec1[0] = (float) (*(poss+(y*xlen+x)*3) +
					gsl_vector_get(vars, 0));
			tvec1[1] = (float) (*(poss+(y*xlen+x)*3+1) +
					gsl_vector_get(vars, 1));
			tvec1[2] = (float) (*(poss+(y*xlen+x)*3+2) +
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


double paraboloid(const double x, const double y, const double f1,
		const double f2) {
	return x*x/(4.0*f1) + y*y/(4.0*f2);
}

double sphere(const double x, const double y, const double rad) {
	// Need to make sure that region dosen't drift outside radius of circle.
	return sqrt(rad*rad - x*x -	y*y);
}

void sphereslope(const float x, const float y, const float rad, float *nrm) {
	// Normal is equal to sphere centre minus sphere surface location.
	*nrm = -x;
	*(nrm+1) = -y;
	*(nrm+2) = -sqrtf(rad*rad - x*x - y*y);
	fscale(1.0f/fnorm(nrm), nrm);
}

void parabslope(const float x, const float y, const float f1, const float f2,
		float *nrm) {
	// Upside down parabaloid, hence the signs.
	*nrm = -x/(2.0f*f1);
	*(nrm+1) = -y/(2.0f*f2);
	*(nrm+2) = -1.0f;
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
			//printf("%10.3e ", gsl_vector_get(minzer->x, i));
		}
		iter++;
	} while (status == GSL_CONTINUE && iter < 2000);

	printf("f() = %7.3f size = %.3f\n", minzer->fval, size);

	for (int i=0; i<n; i++) {
		gsl_vector_set(vars, i, gsl_vector_get(minzer->x, i));
	}

	gsl_multimin_fminimizer_free(minzer);
	minzer = NULL;
}
