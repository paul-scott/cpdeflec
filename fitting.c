#include "fitting.h"

double sqerr(const gsl_vector *vars, void *params) {
	const Errparams *par = (Errparams *) params;
	float *poss = par->poss;
	int xlen = par->xlen;
	int *yb = par->yb;
	double err = 0.0f;

	for (int x=0; x<xlen; x++) {
		for (int y=*(yb+x*2); y<=*(yb+x*2+1); y++) {
			err = err + pow(*(poss+(y*xlen+x)*3+2) -
					(par->f)(*(poss+(y*xlen+x)*3), *(poss+(y*xlen+x)*3+1),
						vars, par->fpars),2.0);
		}
	}
	return err;
}

// This function only works for a fixed curvature paraboloid...
// NEED ROTATIONS 4 THIS.
double fixedparab(const float x, const float y, const gsl_vector *vars,
		const float *as) {
	// Negatives used because needs to be upside down...
	// Could perhaps speed up by removing pow function.
	return -pow((x-gsl_vector_get(vars,0))/(*as),2.0) -
			pow((y-gsl_vector_get(vars,1))/(*(as+1)),2.0) +
			gsl_vector_get(vars,2);
}


double sphere(const float x, const float y, const gsl_vector *vars,
		const float *as) {
	// Need to make sure that region dosen't drift outside radius of circle.
	// z shifted accounting to radius.
	return sqrt(pow(gsl_vector_get(vars,3),2.0) -
			pow(x-gsl_vector_get(vars,0),2.0) -
			pow(y-gsl_vector_get(vars,1),2.0)) - gsl_vector_get(vars,3) +
			gsl_vector_get(vars,2);
}

void minerror(const Errparams *pars, const size_t n, gsl_vector *vars,
		const gsl_vector *stepsize) {
	gsl_multimin_fminimizer *minzer;
	gsl_multimin_function multifunc;
	// Allocate memory for minimizer and set to nelder-mead simplex.
	minzer = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,
			n);
	if (minzer == NULL) {
		printf("Cannot allocate memory for curve fitting...");
		exit(1);
	}

	// Set up multifunc.
	multifunc.n = n;
	multifunc.f = &sqerr;
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

		printf("%5d: ", iter);
		for (int i=0; i<n; i++) {
			printf("%10.3e ", gsl_vector_get(minzer->x, i));
		}
		printf("f() = %7.3f size = %.3f\n", minzer->fval, size);
		iter++;
	} while (status == GSL_CONTINUE && iter < 400);

	for (int i=0; i<n; i++) {
		gsl_vector_set(vars, i, gsl_vector_get(minzer->x, i));
	}

	gsl_multimin_fminimizer_free(minzer);
	minzer = NULL;
}
