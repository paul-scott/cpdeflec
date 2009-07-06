#include "camera.h"

/* Local Structures
 */
typedef struct {
	float *dsep; // Need to change to double??
	float *ph;
} Funcpars;

/* Camera parameters
 * *****************
 * Get the principal values from photogrammetry calibration of camera.
 * Need to sort out coord system of principal point from photogrammetry.
 * For now assume that prpoint is position where optical axis strikes array
 * relative to centre of array. And that increase in direction related to px
 * corresponds in increase in ax dir.
 */
// Can probably pull w and h in from profile.c:
static int camdims[2] = {4288,2848}; // Width and height of CCD array in pixels
static float pxsize = 0.00554f; // Size of a pixel in mm
static float prdist = 20.52f; // Pricipal distance of camera lens in mm
static float prpoint[2] = {0.1065f,-0.2374f}; // Pricipal point in mm

static float campos[3] = {0.0f,0.0f,0.0f};
static float *camtrans;

void print_state(size_t iter, gsl_multiroot_fsolver *s);
void finddirpixcam(float *pix, float *dir);
int func(const gsl_vector *p, void *params, gsl_vector *f);

void initcamera() {
	camtrans = malloc(3*3*sizeof(*camtrans));
}

void freecamera() {
	free(camtrans);
	camtrans = NULL;
}

int objpixsize(float objsize, float objdist) {
	return (int) ceilf((objsize*prdist/objdist)/pxsize);
}

int locatecam(float *dots, float *dotsep, float distguess) {
	/* Assumes first point is origin and second is global y direction from
	 * origin. Global x lies in plane of points 1-2-3.
	 */
	float *gxh = malloc(3*sizeof(*gxh));
	float *gyh = malloc(3*sizeof(*gyh));
	float *gzh = malloc(3*sizeof(*gzh));
	float *fdists = malloc(3*sizeof(*fdists));

	gsl_multiroot_fsolver *solver;
	gsl_multiroot_function multif;
	Funcpars params;
	gsl_vector *dists = gsl_vector_alloc(3);
	float *ph = malloc(3*3*sizeof(*ph));
	int status;
	size_t iter = 0;

	gsl_vector_set_all(dists, (double) distguess);

	// Find directions from pixel values.
	for (int i=0; i<3; i++) finddirpixcam((dots+i*2), (ph+i*3));

	params.dsep = dotsep;
	params.ph = ph;
	// Set up multiroot function.
	multif.f = &func; // System of 3 functions address
	multif.n = 3; // Number of functions
	multif.params = &params; // Parameters to feed func

	// Create solver.
	// Other solver options are _dnewton, _broyden, _hybrid.
	solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, 3);
	if (solver == NULL) {
		printf("Cannot allocate memory for camera locator, exiting...");
		exit(1);
	}
	// Give solver functions and starting vector.
	gsl_multiroot_fsolver_set(solver, &multif, dists);

	do {
		status = gsl_multiroot_fsolver_iterate(solver);

		if (status == GSL_ENOPROG) {
			printf("Solver stuck, exiting...\n");
			exit(1);
		} else if (status == GSL_EBADFUNC) {
			printf("Bad function solver evaluation, exiting...\n");
			exit(1);
		}

		print_state(iter, solver);
 
		status = gsl_multiroot_test_residual(solver->f, 1e-7);

		iter++;
	} while (status == GSL_CONTINUE && iter < 1000);

	*(fdists) = (float) gsl_vector_get(solver->x, 0);
	*(fdists+1) = (float) gsl_vector_get(solver->x, 1);
	*(fdists+2) = (float) gsl_vector_get(solver->x, 2);

	gsl_multiroot_fsolver_free(solver);
	solver = NULL;
	
	// Find unit vectors of global coordinate system described in camera
	// array system.
	// Set dot direction vectors to positions.
	fscale(*(fdists), ph);
	fscale(*(fdists+1), ph+1*3);
	fscale(*(fdists+2), ph+2*3);
	
	// Global y is from dot 0 to dot 1.
	*(gyh) = *(ph+1*3) - *(ph);
	*(gyh+1) = *(ph+1*3+1) - *(ph+1);
	*(gyh+2) = *(ph+1*3+2) - *(ph+2);
	fscale(1.0f/fnorm(gyh), gyh);

	// Appoximate x axis.
	*(gxh) = *(ph+2*3) - *(ph);
	*(gxh+1) = *(ph+2*3+1) - *(ph+1);
	*(gxh+2) = *(ph+2*3+2) - *(ph+2);

	// Find z axis.
	fcross(gxh, gyh, gzh);
	fscale(1.0f/fnorm(gzh), gzh);
	
	// Find x axis.
	fcross(gyh, gzh, gxh);

	// Set camera transform values.
	*(camtrans) = *(gxh);
	*(camtrans+1) = *(gxh+1);
	*(camtrans+2) = *(gxh+2);
	*(camtrans+1*3) = *(gyh);
	*(camtrans+1*3+1) = *(gyh+1);
	*(camtrans+1*3+2) = *(gyh+2);
	*(camtrans+2*3) = *(gzh);
	*(camtrans+2*3+1) = *(gzh+1);
	*(camtrans+2*3+2) = *(gzh+2);

	// Set first ph vector to cam position.
	fscale(-1.0f, ph); // Only uses first vector in ph
	fmatxvec(camtrans, ph, campos); // Only uses first vector in ph

	printf("campos: %f, %f, %f\n", *campos, *(campos+1), *(campos+2));

	gsl_vector_free(dists);
	dists = NULL;

	free(ph);
	ph = NULL;

	free(gxh);
	gxh = NULL;
	free(gyh);
	gyh = NULL;
	free(gzh);
	gzh = NULL;
	free(fdists);
	fdists = NULL;

	return 0; // Return zero if successful.
}

// Finds pixel that intercepts position vector.
void findpix(float *vec, int *pix) {

}

// Calculates global direction of pixel from camera.
void finddirpix(int *pix, float *dir) {

	// Note that finddirpixcam requires a float...
	// Might consider making another function that only operates on integers.
}

// Calculates camera coord direction of pixel from camera. Note the use of
// float to get sub-pixel accuracy for locating.
void finddirpixcam(float *pix, float *dir) {
	/* Define camera ax opposite to px and ay opposite to py. az is along
	 * optical axis, therefore (x,y,z)=(0,0,0) optical centre. As position
	 * of pixel in array increases, image pixel value increases.
	 * Signs used here to get the direction right. Reverse signs and we
	 * get position of pixel in the array from optical centre.
	 */
	float mag;
	*(dir) = -(pxsize*(*(pix) - *(camdims)/2.0f + 0.5f) + *(prpoint));
	*(dir+1) = -(pxsize*(*(pix+1) - *(camdims+1)/2.0f + 0.5f) + *(prpoint+1));
	*(dir+2) = prdist;
	// Normalizing vector.
	fscale(1.0f/fnorm(dir),dir);
}

int func(const gsl_vector *p, void *params, gsl_vector *f) {
	Funcpars *par = (Funcpars *) params;
	const double p0 = gsl_vector_get(p,0);
	const double p1 = gsl_vector_get(p,1);
	const double p2 = gsl_vector_get(p,2);
	// NEED TO CHECK WHAT HAPPENS WHEN WE MULTIPLY A FLOAT WITH DOUBLE.
	double f1 = p0*p0-2*p0*p1*fdot((par->ph),(par->ph+1*3))+p1*p1-
			(*(par->dsep))*(*(par->dsep));
	double f2 = p1*p1-2*p1*p2*fdot((par->ph+1*3),(par->ph+2*3))+p2*p2-
			(*(par->dsep+1))*(*(par->dsep+1));
	double f3 = p2*p2-2*p2*p0*fdot((par->ph+2*3),(par->ph))+p0*p0-
			(*(par->dsep+2))*(*(par->dsep+2));

	gsl_vector_set(f, 0, f1);
	gsl_vector_set(f, 1, f2);
	gsl_vector_set(f, 2, f3);

	return GSL_SUCCESS;
}

void print_state(size_t iter, gsl_multiroot_fsolver *s) {
	// Prints current state of solver.
	printf("iter = %3u x = %.3f %.3f %.3f f(x) = %.3e %.3e %.3e\n",
			iter,
			gsl_vector_get(s->x, 0),
			gsl_vector_get(s->x, 1),
			gsl_vector_get(s->x, 2),
			gsl_vector_get(s->f, 0),
			gsl_vector_get(s->f, 1),
			gsl_vector_get(s->f, 2));
}

